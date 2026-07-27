#include "mex.h"
#include "matrix.h"
#include "lapack.h"

#include <math.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>

/*
 * Persistent complex Hermitian Cholesky storage in LAPACK RFP format.
 *
 * Compile with the interleaved-complex API and MATLAB LAPACK/BLAS:
 *   mex -R2018a rfp_chol_mex.c -lmwlapack -lmwblas
 *
 * Commands:
 *   rfp_chol_mex('init', n)
 *   rfp_chol_mex('put', rows, cols, block)
 *   rfp_chol_mex('factor')
 *   X = rfp_chol_mex('solve', B)
 *   s = rfp_chol_mex('info')
 *   rfp_chol_mex('free')
 *
 * rows and cols are one-based MATLAB indices.  put stores only coordinates
 * satisfying rows(i) <= cols(j); entries below the diagonal in block are
 * intentionally ignored.  The caller must supply the final upper-Hermitian
 * values.  Repeated coordinates overwrite earlier values.
 */

typedef char complex_size_check[
    (sizeof(mxComplexDouble) == 2U * sizeof(double)) ? 1 : -1];
typedef char complex_imag_offset_check[
    (offsetof(mxComplexDouble, imag) == sizeof(double)) ? 1 : -1];

static mxComplexDouble *g_arf = NULL;
static mwSize g_n = 0;
static mwSize g_nt = 0;
static uint64_T g_put_calls = 0;
static uint64_T g_written_entries = 0;
static uint64_T g_ignored_lower_entries = 0;
static ptrdiff_t g_lapack_info = 0;
static int g_factor_attempted = 0;
static int g_factorized = 0;
static int g_locked = 0;
static int g_exit_registered = 0;

/* Reset the persistent matrix metadata. */
static void reset_state(void)
{
    g_n = 0;
    g_nt = 0;
    g_put_calls = 0;
    g_written_entries = 0;
    g_ignored_lower_entries = 0;
    g_lapack_info = 0;
    g_factor_attempted = 0;
    g_factorized = 0;
}

/* Release the packed matrix buffer and its metadata. */
static void release_buffer(void)
{
    if (g_arf != NULL) {
        mxFree(g_arf);
        g_arf = NULL;
    }
    reset_state();
}

/* Release persistent state when MATLAB unloads the MEX file. */
static void cleanup_at_exit(void)
{
    release_buffer();
    g_locked = 0;
}

/* Verify that packed matrix storage has been initialized. */
static void require_initialized(void)
{
    if (g_arf == NULL || g_n == 0) {
        mexErrMsgIdAndTxt("Example4:hartreeComparison:rfp:notInitialized",
                          "Call rfp_chol_mex('init', n) first.");
    }
}

/* Validate the number of command inputs. */
static void require_nrhs(int nrhs, int expected, const char *command)
{
    if (nrhs != expected) {
        mexErrMsgIdAndTxt("Example4:hartreeComparison:rfp:inputCount",
                          "Command '%s' requires %d input argument(s) after the command.",
                          command, expected - 1);
    }
}

/* Validate the maximum number of command outputs. */
static void require_nlhs_at_most(int nlhs, int maximum, const char *command)
{
    if (nlhs > maximum) {
        mexErrMsgIdAndTxt("Example4:hartreeComparison:rfp:outputCount",
                          "Command '%s' returns at most %d output argument(s).",
                          command, maximum);
    }
}

/* Convert the matrix order to a checked MATLAB size. */
static mwSize checked_order(const mxArray *value)
{
    double n_value;
    mwSize n;

    if (!mxIsDouble(value) || mxIsComplex(value) || mxIsSparse(value) ||
        mxGetNumberOfElements(value) != 1) {
        mexErrMsgIdAndTxt("Example4:hartreeComparison:rfp:invalidOrder",
                          "n must be a real, full, double scalar.");
    }

    n_value = mxGetScalar(value);
    if (!mxIsFinite(n_value) || n_value < 1.0 || floor(n_value) != n_value ||
        n_value > (double)PTRDIFF_MAX || n_value > (double)SIZE_MAX) {
        mexErrMsgIdAndTxt("Example4:hartreeComparison:rfp:invalidOrder",
                          "n must be a positive integer representable by LAPACK.");
    }

    n = (mwSize)n_value;
    if ((double)n != n_value) {
        mexErrMsgIdAndTxt("Example4:hartreeComparison:rfp:invalidOrder",
                          "n cannot be represented exactly by mwSize.");
    }
    return n;
}

/* Compute the number of entries in packed triangular storage. */
static mwSize packed_length(mwSize n)
{
    mwSize a;
    mwSize b;

    if ((n & 1U) == 0U) {
        a = n / 2U;
        b = n + 1U;
    } else {
        a = n;
        b = (n + 1U) / 2U;
    }

    if (a != 0U && b > SIZE_MAX / a) {
        mexErrMsgIdAndTxt("Example4:hartreeComparison:rfp:sizeOverflow",
                          "n*(n+1)/2 overflows the addressable size.");
    }
    return a * b;
}

/* Validate one-based row or column index vectors. */
static void validate_index_vector(const mxArray *array, const char *name)
{
    const double *values;
    mwSize count;
    mwSize index;

    if (!mxIsDouble(array) || mxIsComplex(array) || mxIsSparse(array) ||
        mxGetNumberOfDimensions(array) > 2) {
        mexErrMsgIdAndTxt("Example4:hartreeComparison:rfp:invalidIndices",
                          "%s must be a real, full, double vector.", name);
    }

    count = mxGetNumberOfElements(array);
    if (count != 0U && mxGetM(array) != 1U && mxGetN(array) != 1U) {
        mexErrMsgIdAndTxt("Example4:hartreeComparison:rfp:invalidIndices",
                          "%s must be a vector.", name);
    }

    values = mxGetDoubles(array);
    for (index = 0; index < count; ++index) {
        double value = values[index];
        if (!mxIsFinite(value) || value < 1.0 || value > (double)g_n ||
            floor(value) != value || (double)(mwIndex)value != value) {
            mexErrMsgIdAndTxt("Example4:hartreeComparison:rfp:indexOutOfRange",
                              "%s contains an invalid matrix index.", name);
        }
    }
}

/*
 * Map zero-based A(row,col), row <= col, to TRANSR='N', UPLO='U' RFP.
 * The first floor(n/2) columns of upper A occupy the conjugated triangle
 * below the rectangular trapezoid, hence conjugate_store is set there.
 */
static mwIndex upper_rfp_offset(mwIndex row, mwIndex col,
                                int *conjugate_store)
{
    mwIndex k = (mwIndex)(g_n / 2U);
    mwIndex leading_dimension = (mwIndex)(g_n + ((g_n & 1U) == 0U));
    mwIndex rfp_row;
    mwIndex rfp_col;

    if (col >= k) {
        rfp_row = row;
        rfp_col = col - k;
        *conjugate_store = 0;
    } else {
        rfp_row = k + 1U + col;
        rfp_col = row;
        *conjugate_store = 1;
    }

    return rfp_row + leading_dimension * rfp_col;
}

/* Package the current persistent-state statistics. */
static mxArray *make_info_struct(void)
{
    const char *fields[] = {
        "initialized", "n", "elements", "bytes", "factorAttempted",
        "factorized", "lapackInfo", "putCalls", "writtenEntries",
        "ignoredLowerEntries", "transr", "uplo"
    };
    mxArray *result = mxCreateStructMatrix(1, 1, 12, fields);

    mxSetField(result, 0, "initialized", mxCreateLogicalScalar(g_arf != NULL));
    mxSetField(result, 0, "n", mxCreateDoubleScalar((double)g_n));
    mxSetField(result, 0, "elements", mxCreateDoubleScalar((double)g_nt));
    mxSetField(result, 0, "bytes",
               mxCreateDoubleScalar((double)g_nt * (double)sizeof(mxComplexDouble)));
    mxSetField(result, 0, "factorAttempted",
               mxCreateLogicalScalar(g_factor_attempted != 0));
    mxSetField(result, 0, "factorized", mxCreateLogicalScalar(g_factorized != 0));
    mxSetField(result, 0, "lapackInfo", mxCreateDoubleScalar((double)g_lapack_info));
    mxSetField(result, 0, "putCalls", mxCreateDoubleScalar((double)g_put_calls));
    mxSetField(result, 0, "writtenEntries",
               mxCreateDoubleScalar((double)g_written_entries));
    mxSetField(result, 0, "ignoredLowerEntries",
               mxCreateDoubleScalar((double)g_ignored_lower_entries));
    mxSetField(result, 0, "transr", mxCreateString("N"));
    mxSetField(result, 0, "uplo", mxCreateString("U"));
    return result;
}

/* Allocate persistent RFP storage for a matrix order. */
static void command_init(int nlhs, mxArray *plhs[], int nrhs,
                         const mxArray *prhs[])
{
    mwSize n;
    mwSize nt;

    require_nrhs(nrhs, 2, "init");
    require_nlhs_at_most(nlhs, 1, "init");
    if (g_arf != NULL) {
        mexErrMsgIdAndTxt("Example4:hartreeComparison:rfp:alreadyInitialized",
                          "The RFP buffer already exists; call 'free' first.");
    }

    n = checked_order(prhs[1]);
    nt = packed_length(n);
    if (nt > SIZE_MAX / sizeof(mxComplexDouble)) {
        mexErrMsgIdAndTxt("Example4:hartreeComparison:rfp:sizeOverflow",
                          "The complex RFP buffer size overflows size_t.");
    }

    g_arf = (mxComplexDouble *)mxCalloc(nt, sizeof(mxComplexDouble));
    if (g_arf == NULL) {
        mexErrMsgIdAndTxt("Example4:hartreeComparison:rfp:allocationFailed",
                          "MATLAB could not allocate the complex RFP buffer.");
    }
    mexMakeMemoryPersistent(g_arf);
    g_n = n;
    g_nt = nt;
    g_lapack_info = 0;
    g_factor_attempted = 0;
    g_factorized = 0;
    g_put_calls = 0;
    g_written_entries = 0;
    g_ignored_lower_entries = 0;

    if (!g_locked) {
        mexLock();
        g_locked = 1;
    }

    if (nlhs == 1) {
        plhs[0] = make_info_struct();
    }
}

/* Copy an upper-triangular block into RFP storage. */
static void command_put(int nlhs, int nrhs, const mxArray *prhs[])
{
    const mxArray *rows_array;
    const mxArray *cols_array;
    const mxArray *block_array;
    const double *rows;
    const double *cols;
    const double *real_block = NULL;
    const mxComplexDouble *complex_block = NULL;
    mwSize row_count;
    mwSize col_count;
    mwIndex block_col;

    require_nrhs(nrhs, 4, "put");
    require_nlhs_at_most(nlhs, 0, "put");
    require_initialized();
    if (g_factor_attempted) {
        mexErrMsgIdAndTxt("Example4:hartreeComparison:rfp:factorAlreadyAttempted",
                          "The buffer cannot be modified after factor() was attempted.");
    }

    rows_array = prhs[1];
    cols_array = prhs[2];
    block_array = prhs[3];
    validate_index_vector(rows_array, "rows");
    validate_index_vector(cols_array, "cols");
    row_count = mxGetNumberOfElements(rows_array);
    col_count = mxGetNumberOfElements(cols_array);

    if (!mxIsDouble(block_array) || mxIsSparse(block_array) ||
        mxGetNumberOfDimensions(block_array) > 2 ||
        mxGetM(block_array) != row_count || mxGetN(block_array) != col_count) {
        mexErrMsgIdAndTxt("Example4:hartreeComparison:rfp:invalidBlock",
                          "block must be a full double matrix of size numel(rows)-by-numel(cols).");
    }

    rows = mxGetDoubles(rows_array);
    cols = mxGetDoubles(cols_array);
    if (mxIsComplex(block_array)) {
        complex_block = mxGetComplexDoubles(block_array);
    } else {
        real_block = mxGetDoubles(block_array);
    }

    for (block_col = 0; block_col < col_count; ++block_col) {
        mwIndex col = (mwIndex)cols[block_col] - 1U;
        mwIndex block_row;
        for (block_row = 0; block_row < row_count; ++block_row) {
            mwIndex row = (mwIndex)rows[block_row] - 1U;
            mwIndex source_offset = block_row + row_count * block_col;
            mxComplexDouble value;
            mwIndex destination_offset;
            int conjugate_store;

            if (row > col) {
                ++g_ignored_lower_entries;
                continue;
            }

            if (complex_block != NULL) {
                value = complex_block[source_offset];
            } else {
                value.real = real_block[source_offset];
                value.imag = 0.0;
            }

            if (row == col) {
                value.imag = 0.0;
            }

            destination_offset = upper_rfp_offset(row, col, &conjugate_store);
            if (destination_offset >= g_nt) {
                mexErrMsgIdAndTxt("Example4:hartreeComparison:rfp:internalMapping",
                                  "Internal RFP coordinate mapping exceeded the buffer.");
            }
            if (conjugate_store) {
                value.imag = -value.imag;
            }
            g_arf[destination_offset] = value;
            ++g_written_entries;
        }
    }
    ++g_put_calls;
}

/* Factor the packed Hermitian matrix with LAPACK. */
static void command_factor(int nlhs, mxArray *plhs[], int nrhs)
{
    const char transr = 'N';
    const char uplo = 'U';
    ptrdiff_t n_lapack;
    ptrdiff_t info = 0;

    require_nrhs(nrhs, 1, "factor");
    require_nlhs_at_most(nlhs, 1, "factor");
    require_initialized();
    if (g_factor_attempted) {
        mexErrMsgIdAndTxt("Example4:hartreeComparison:rfp:factorAlreadyAttempted",
                          "factor() may be called only once for each init().");
    }

    n_lapack = (ptrdiff_t)g_n;
    g_factor_attempted = 1;
    zpftrf(&transr, &uplo, &n_lapack, (double *)(void *)g_arf, &info);
    g_lapack_info = info;
    g_factorized = (info == 0);

    if (nlhs == 1) {
        plhs[0] = mxCreateDoubleScalar((double)info);
    }
    if (info < 0) {
        mexErrMsgIdAndTxt("Example4:hartreeComparison:rfp:lapackArgument",
                          "zpftrf rejected LAPACK argument %td.", -info);
    }
    if (info > 0) {
        mexErrMsgIdAndTxt("Example4:hartreeComparison:rfp:notPositiveDefinite",
                          "zpftrf found a non-positive leading minor at index %td.", info);
    }
}

/* Solve the factorized system for one or more right-hand sides. */
static void command_solve(int nlhs, mxArray *plhs[], int nrhs,
                          const mxArray *prhs[])
{
    const mxArray *right_hand_side;
    const double *real_input = NULL;
    const mxComplexDouble *complex_input = NULL;
    mxArray *solution;
    mxComplexDouble *solution_values;
    mwSize element_count;
    mwSize index;
    mwSize rhs_count;
    ptrdiff_t n_lapack;
    ptrdiff_t nrhs_lapack;
    ptrdiff_t ldb;
    ptrdiff_t info = 0;
    const char transr = 'N';
    const char uplo = 'U';

    require_nrhs(nrhs, 2, "solve");
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("Example4:hartreeComparison:rfp:outputCount",
                          "solve requires exactly one output argument.");
    }
    require_initialized();
    if (!g_factorized) {
        mexErrMsgIdAndTxt("Example4:hartreeComparison:rfp:notFactorized",
                          "A successful factor() call is required before solve().");
    }

    right_hand_side = prhs[1];
    if (!mxIsDouble(right_hand_side) || mxIsSparse(right_hand_side) ||
        mxGetNumberOfDimensions(right_hand_side) > 2 ||
        mxGetM(right_hand_side) != g_n) {
        mexErrMsgIdAndTxt("Example4:hartreeComparison:rfp:invalidRightHandSide",
                          "B must be a full double matrix with n rows.");
    }

    rhs_count = mxGetN(right_hand_side);
    if (rhs_count > (mwSize)PTRDIFF_MAX) {
        mexErrMsgIdAndTxt("Example4:hartreeComparison:rfp:tooManyRightHandSides",
                          "The number of right-hand sides exceeds LAPACK integer range.");
    }
    if (rhs_count != 0U && g_n > SIZE_MAX / rhs_count) {
        mexErrMsgIdAndTxt("Example4:hartreeComparison:rfp:rightHandSideOverflow",
                          "The right-hand-side element count overflows size_t.");
    }
    element_count = g_n * rhs_count;
    solution = mxCreateDoubleMatrix(g_n, rhs_count, mxCOMPLEX);
    solution_values = mxGetComplexDoubles(solution);

    if (mxIsComplex(right_hand_side)) {
        complex_input = mxGetComplexDoubles(right_hand_side);
        memcpy(solution_values, complex_input,
               element_count * sizeof(mxComplexDouble));
    } else {
        real_input = mxGetDoubles(right_hand_side);
        for (index = 0; index < element_count; ++index) {
            solution_values[index].real = real_input[index];
            solution_values[index].imag = 0.0;
        }
    }

    n_lapack = (ptrdiff_t)g_n;
    nrhs_lapack = (ptrdiff_t)rhs_count;
    ldb = n_lapack;
    zpftrs(&transr, &uplo, &n_lapack, &nrhs_lapack,
            (const double *)(const void *)g_arf,
            (double *)(void *)solution_values, &ldb, &info);
    if (info != 0) {
        mxDestroyArray(solution);
        mexErrMsgIdAndTxt("Example4:hartreeComparison:rfp:lapackSolve",
                          "zpftrs returned LAPACK info %td.", info);
    }
    plhs[0] = solution;
}

/* Release the persistent RFP matrix. */
static void command_free(int nlhs, int nrhs)
{
    require_nrhs(nrhs, 1, "free");
    require_nlhs_at_most(nlhs, 0, "free");
    release_buffer();
    if (g_locked) {
        g_locked = 0;
        mexUnlock();
    }
}

/* Return the persistent matrix state and counters. */
static void command_info(int nlhs, mxArray *plhs[], int nrhs)
{
    require_nrhs(nrhs, 1, "info");
    require_nlhs_at_most(nlhs, 1, "info");
    if (nlhs == 1) {
        plhs[0] = make_info_struct();
    }
}

/* Dispatch MATLAB commands to the packed-matrix operations. */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs,
                 const mxArray *prhs[])
{
    char command[32];

    if (!g_exit_registered) {
        mexAtExit(cleanup_at_exit);
        g_exit_registered = 1;
    }

    if (nrhs < 1 || !mxIsChar(prhs[0]) ||
        mxGetString(prhs[0], command, sizeof(command)) != 0) {
        mexErrMsgIdAndTxt("Example4:hartreeComparison:rfp:invalidCommand",
                          "The first input must be a command character vector.");
    }

    if (strcmp(command, "init") == 0) {
        command_init(nlhs, plhs, nrhs, prhs);
    } else if (strcmp(command, "put") == 0) {
        command_put(nlhs, nrhs, prhs);
    } else if (strcmp(command, "factor") == 0) {
        command_factor(nlhs, plhs, nrhs);
    } else if (strcmp(command, "solve") == 0) {
        command_solve(nlhs, plhs, nrhs, prhs);
    } else if (strcmp(command, "free") == 0) {
        command_free(nlhs, nrhs);
    } else if (strcmp(command, "info") == 0) {
        command_info(nlhs, plhs, nrhs);
    } else {
        mexErrMsgIdAndTxt("Example4:hartreeComparison:rfp:unknownCommand",
                          "Unknown command '%s'.", command);
    }
}
