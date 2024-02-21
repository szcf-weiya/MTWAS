#include<Rcpp.h>
const double pi = 3.14159265358979323846;
using namespace Rcpp;
#include<RcppGSL.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_combination.h>
#include<gsl/gsl_statistics.h>
#include<gsl/gsl_fit.h>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_multifit.h>
using namespace std;

double BIC(int N, int p, double rss) {
    double ll = -0.5 * N * (log(2*pi) + 1 - log(N) + log(rss));
    return -2 * ll + log(N) * (p + 1);
}

// adapted from "gsl-2.7/multifit/linear_common.c"
static int
    multifit_linear_solve (const gsl_matrix * X,
                           const gsl_vector * y,
                           const double tol,
                           const double lambda,
                           size_t * rank,
                           gsl_vector * c,
                           double *rnorm,
                           double *snorm,
                           gsl_multifit_linear_workspace * work)
    {
        const size_t n = X->size1;
        const size_t p = X->size2;

        if (n != work->n || p != work->p)
        {
            GSL_ERROR("observation matrix does not match workspace", GSL_EBADLEN);
        }
        else if (n != y->size)
        {
            GSL_ERROR("number of observations in y does not match matrix",
                      GSL_EBADLEN);
        }
        else if (p != c->size)
        {
            GSL_ERROR ("number of parameters c does not match matrix",
                       GSL_EBADLEN);
        }
        else if (tol <= 0)
        {
            GSL_ERROR ("tolerance must be positive", GSL_EINVAL);
        }
        else
        {
            const double lambda_sq = lambda * lambda;

            double rho_ls = 0.0;     /* contribution to rnorm from OLS */

            size_t j, p_eff;

            /* these inputs are previously computed by multifit_linear_svd() */
            gsl_matrix_view A = gsl_matrix_submatrix(work->A, 0, 0, n, p);
            gsl_matrix_view Q = gsl_matrix_submatrix(work->Q, 0, 0, p, p);
            gsl_vector_view S = gsl_vector_subvector(work->S, 0, p);

            /* workspace */
            gsl_matrix_view QSI = gsl_matrix_submatrix(work->QSI, 0, 0, p, p);
            gsl_vector_view xt = gsl_vector_subvector(work->xt, 0, p);
            gsl_vector_view D = gsl_vector_subvector(work->D, 0, p);
            gsl_vector_view t = gsl_vector_subvector(work->t, 0, n);

            /*
             * Solve y = A c for c
             * c = Q diag(s_i / (s_i^2 + lambda_i^2)) U^T y
             */

            /* compute xt = U^T y */
            gsl_blas_dgemv (CblasTrans, 1.0, &A.matrix, y, 0.0, &xt.vector);

            if (n > p)
            {
                /*
                 * compute OLS residual norm = || y - U U^T y ||;
                 * for n = p, U U^T = I, so no need to calculate norm
                 */
                gsl_vector_memcpy(&t.vector, y);
                gsl_blas_dgemv(CblasNoTrans, -1.0, &A.matrix, &xt.vector, 1.0, &t.vector);
                rho_ls = gsl_blas_dnrm2(&t.vector);
            }

            if (lambda > 0.0)
            {
                /* xt <-- [ s(i) / (s(i)^2 + lambda^2) ] .* U^T y */
                for (j = 0; j < p; ++j)
                {
                    double sj = gsl_vector_get(&S.vector, j);
                    double f = (sj * sj) / (sj * sj + lambda_sq);
                    double *ptr = gsl_vector_ptr(&xt.vector, j);

                    /* use D as workspace for residual norm */
                    gsl_vector_set(&D.vector, j, (1.0 - f) * (*ptr));

                    *ptr *= sj / (sj*sj + lambda_sq);
                }

                /* compute regularized solution vector */
                gsl_blas_dgemv (CblasNoTrans, 1.0, &Q.matrix, &xt.vector, 0.0, c);

                /* compute solution norm */
                *snorm = gsl_blas_dnrm2(c);

                /* compute residual norm */
                *rnorm = gsl_blas_dnrm2(&D.vector);

                if (n > p)
                {
                    /* add correction to residual norm (see eqs 6-7 of [1]) */
                    *rnorm = sqrt((*rnorm) * (*rnorm) + rho_ls * rho_ls);
                }

                /* reset D vector */
                gsl_vector_set_all(&D.vector, 1.0);
            }
            else
            {
                /* Scale the matrix Q, QSI = Q S^{-1} */

                gsl_matrix_memcpy (&QSI.matrix, &Q.matrix);

                {
                    double s0 = gsl_vector_get (&S.vector, 0);
                    p_eff = 0;

                    for (j = 0; j < p; j++)
                    {
                        gsl_vector_view column = gsl_matrix_column (&QSI.matrix, j);
                        double sj = gsl_vector_get (&S.vector, j);
                        double alpha;

                        if (sj <= tol * s0)
                        {
                            alpha = 0.0;
                        }
                        else
                        {
                            alpha = 1.0 / sj;
                            p_eff++;
                        }

                        gsl_vector_scale (&column.vector, alpha);
                    }

                    *rank = p_eff;
                }

                gsl_blas_dgemv (CblasNoTrans, 1.0, &QSI.matrix, &xt.vector, 0.0, c);

                /* Unscale the balancing factors */
                gsl_vector_div (c, &D.vector);

                *snorm = gsl_blas_dnrm2(c);

                /// begin of add by weiya [2022-08-29]
                gsl_vector_memcpy(&t.vector, y);
                gsl_blas_dgemv(CblasNoTrans, -1.0, X, c, 1.0, &t.vector);
                rho_ls = gsl_blas_dnrm2(&t.vector); // end of add
                *rnorm = rho_ls;
            }

            return GSL_SUCCESS;
        }
    }

double lm(const RcppGSL::Matrix &X, const RcppGSL::Vector &y) {

    int n = X.nrow(), k = X.ncol();
    double chisq;

    RcppGSL::Vector coef(k);                // to hold the coefficient vector
    RcppGSL::Matrix cov(k,k);               // and the covariance matrix

    // the actual fit requires working memory we allocate and free
    gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc (n, k);
    int status;
    size_t rank;
    double rnorm = 0.0, snorm;
    status = gsl_multifit_linear_bsvd(X, work);
    status = multifit_linear_solve(X, y, GSL_DBL_EPSILON, -1.0, &rank, coef, &rnorm, &snorm, work);
    chisq = rnorm * rnorm;
    // gsl_multifit_linear (X, y, coef, cov, &chisq, work);
    // gsl_multifit_linear_tsvd(X, y, 1e-7, coef, cov, &chisq, &rank, work);
    gsl_multifit_linear_free (work);
    // cout << "rank = " << rank << endl;
    // for (size_t i = 0; i < k; i++)
    //     cout << coef[i] << " ";
    // cout << endl;
    return chisq;
}

double calc_BIC(const RcppGSL::Matrix &xx, const RcppGSL::Vector &yy, const std::vector<int> terms, double gam = 0) {
    int d = terms.size();
    int N = xx.nrow(), D = xx.ncol();
    double rss, bic;
    if (d == 0) {
        rss = gsl_stats_tss(yy->data, 1, N);
        bic = BIC(N, 1, rss) + 2 * gam * log(D);
    } else {
        gsl_matrix *X;
        gsl_vector *tmp;
        X = gsl_matrix_alloc(N, d + 1);
        tmp = gsl_vector_alloc(N);
        for(size_t i = 0; i < d; i++) {
            gsl_matrix_get_col(tmp, xx, terms[i]-1);
            gsl_matrix_set_col(X, i+1, tmp);
        }
        gsl_vector_set_all(tmp, 1);
        gsl_matrix_set_col(X, 0, tmp);
        rss = lm(X, yy);
        bic = BIC(N, d+1, rss) + (1+d) * 2 * gam * log(D);
    }
    return bic;
}

//' Select main effect via SODA
//'
//' @param xx design matrix
//' @param yy response vector
//' @param xnames colnames for the design matrix
//' @param fixset an index vector containing effects are fixed
//' @param norm deprecated (TO REMOVE)
//' @param verbose whether to print info
//' @param gam tuning parameter \eqn{\gamma} in the EBIC criterion
//' @param minF0 the smallest number of main effects
//' @param allow_empty whether to allow empty main effects
//' @return the best EBIC and the corresponding selected set of main effects
//' @export
// [[Rcpp::depends(RcppGSL)]]
// [[Rcpp::export]]
Rcpp::List soda_main_allback(const RcppGSL::Matrix &xx, const RcppGSL::Vector &yy,
                             const CharacterVector &xnames,
                             const std::vector<int> &fixset,
                             bool norm = false, bool verbose = false,
                             double gam = 0, int minF0 = 3,
                             bool allow_empty = false) {
    int N = xx.nrow(), D = xx.ncol();
    minF0 = min(D, minF0);
    vector<double> BICs;
    vector<int> terms = fixset;
    BICs.push_back(calc_BIC(xx, yy, terms, gam));
    if (verbose) Rcout << "Initialization: EBIC = " << BICs[0] << endl;

    vector<int> set_all;
    vector<int> cur_set = fixset;
    for (size_t i = 0; i < D; i++) { set_all.push_back(i+1); }
    double cur_score, new_score;
    if (minF0 == 0)
        cur_score = BICs[1];
    else
        cur_score = 1e6;
    /************************
     * Linear Forward Stage *
     ************************/
    if (verbose) Rcout << "Forward Stage - Main effects:" << endl;
    // int niter = 0;
    while (true) {
        /**********************
         * Forward Operations *
         **********************/
        vector<int> not_set;
        // require set_all and cur_set are sorted
        set_difference(set_all.begin(), set_all.end(), cur_set.begin(), cur_set.end(),
                                                inserter(not_set, not_set.begin()));
        int Nnset = not_set.size();
        if (Nnset == 0 || cur_set.size() >= ((yy.size())/5))  // n_ops == 0
            break;
        vector<int> new_set = cur_set; // does not change cur_set (https://cplusplus.com/reference/vector/vector/operator=/)
        new_set.push_back(-1);
        int ncur = cur_set.size();
        if (ncur < minF0)
            cur_score = 1e6;
        bool flag = false; // there exists smaller new scores;
        double best_score = 1e6 + 1;
        int best_term;
        for (size_t j = 0; j < Nnset; j++) {
            int jj = not_set[j];
            new_set[ncur] = jj;
            vector<int> sort_new_set = new_set;
            // partial_sort_copy(begin(new_set), end(new_set), begin(sort_new_set), end(sort_new_set));
            // sort(sort_new_set.begin(), sort_new_set.end());
            new_score = calc_BIC(xx, yy, new_set, gam);
            if ( new_score < cur_score || ncur < minF0) {
                flag = true;

                /* pick the best */
                if (new_score < best_score) {
                    best_score = new_score;
                    best_term = jj;
                }
            }
        }
        if (!flag)
            break;
        BICs.push_back(best_score);
        cur_set.push_back(best_term);
        cur_score = best_score;
        sort(cur_set.begin(), cur_set.end());
        if (verbose)
            Rcout << "  Main effects: add variable " << best_term << ": " << xnames[best_term-1] << " into selection set... df = " << cur_set.size() + 1 << ", EBIC = " << BICs.back() << endl;
    }
    if (verbose) Rcout << "Backward stage:" << endl;
    /******************
     * Backward Stage *
     ******************/
    // since only linear terms, cur_set is equivalent to cur_terms
    if (cur_set.size() > 12) {
        while(true) {
            /***********************
             * Backward Operations *
             ***********************/
            int Nterms = cur_set.size();
            vector<int> cur_terms_nfix;
            set_difference(cur_set.begin(), cur_set.end(), fixset.begin(), fixset.end(),
                           inserter(cur_terms_nfix, cur_terms_nfix.begin()));
            int Nterms_nfix = cur_terms_nfix.size();
            double best_score = 1e6+1;
            int best_term;
            bool flag = false;
            if (Nterms_nfix > 12) {
                vector<int> new_terms = cur_terms_nfix;
                new_terms.pop_back(); // delete the last element
                // put the fixset at the end
                for(size_t i = 0; i < fixset.size(); i++) {
                    new_terms.push_back(fixset[i]);
                }
                for(size_t j = 0; j < Nterms_nfix; j++){
                    int term = cur_terms_nfix[j];
                    int backup = new_terms[j];
                    if (j != Nterms_nfix - 1)
                        new_terms[j] = cur_terms_nfix.back(); // delete the j-th element by swapping the deleted last element
                    new_score = calc_BIC(xx, yy, new_terms, gam);
                    if ( new_score < cur_score) {
                        flag = true;
                        /* pick the best */
                        if (new_score < best_score) {
                            best_score = new_score;
                            best_term = j;
                        }
                    }
                    if (j != Nterms_nfix - 1) // reset new_terms
                        new_terms[j] = backup;
                }
            }
            if (!flag)
                break;
            BICs.push_back(best_score);
            int best_term_idx = cur_terms_nfix[best_term];
            cur_set= cur_terms_nfix;
            cur_set.erase(cur_set.begin()+best_term);
            for (size_t i = 0; i < fixset.size(); i++)
                cur_set.push_back(fixset[i]);
            sort(cur_set.begin(), cur_set.end());
            cur_score = best_score;
            if (verbose)
                Rcout << "  Remove term " << xnames[best_term_idx-1] << ": " << "from selection set...  df = " << cur_set.size() + 1 << ", EBIC = " << BICs.back() << endl;
        }
    }
    gsl_combination *c;
    vector<int> cur_terms_nfix;
    set_difference(cur_set.begin(), cur_set.end(), fixset.begin(), fixset.end(),
                   inserter(cur_terms_nfix, cur_terms_nfix.begin()));
    int ncur = cur_terms_nfix.size();
    double best_score = 1e7;
    vector<int> best_set;
    double tmp_score;
    for (int i = (int)!allow_empty; i <= ncur; i++) {
        c = gsl_combination_calloc(ncur, i);
        do {
            vector<int> term;
            for (size_t j = 0; j < i; j++) {
                term.push_back(cur_terms_nfix[gsl_combination_get(c, j)]);
            }
            for (size_t j = 0; j < fixset.size(); j++)
                term.push_back(fixset[j]);
            tmp_score = calc_BIC(xx, yy, term, gam);
            if (tmp_score < best_score) {
                best_score = tmp_score;
                best_set = term;
            }
        } while (gsl_combination_next(c) == GSL_SUCCESS);
    }
    gsl_combination_free(c);
    if (verbose) {
        Rcout << "       terms: ";
        for (size_t i = 0; i < best_set.size(); i++)
            Rcout << xnames[best_set[i]-1] << " ";
        Rcout << endl;
    }
    return Rcpp::List::create(Rcpp::Named("EBIC") = best_score,
                              Rcpp::Named("Select Set") = best_set);
}

/*** R
# system.time(
# for (i in 1:20)
#     #print(soda_main_allback(xx0[[i]], yy0[[i]], colnames(xx0[[i]]), FALSE, FALSE, 1, 3))
#     soda_main_allback(xx0[[i]], yy0[[i]], colnames(xx0[[i]]), FALSE, FALSE, 1, 3)
# )
# soda_main_allback(train.x, train.y, 1:ncol(train.x), numeric(0), FALSE, FALSE, 1, 16)
*/
