// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// find_seqnum
IntegerVector find_seqnum(IntegerVector clens, IntegerVector idx);
RcppExport SEXP _neuroim2_find_seqnum(SEXP clensSEXP, SEXP idxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type clens(clensSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type idx(idxSEXP);
    rcpp_result_gen = Rcpp::wrap(find_seqnum(clens, idx));
    return rcpp_result_gen;
END_RCPP
}
// grid_to_intvec
int grid_to_intvec(IntegerVector D, IntegerVector vox);
RcppExport SEXP _neuroim2_grid_to_intvec(SEXP DSEXP, SEXP voxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type D(DSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type vox(voxSEXP);
    rcpp_result_gen = Rcpp::wrap(grid_to_intvec(D, vox));
    return rcpp_result_gen;
END_RCPP
}
// exgridToIndex3DCpp
IntegerVector exgridToIndex3DCpp(IntegerVector array_dim, IntegerVector iind, IntegerVector jind, IntegerVector kind);
RcppExport SEXP _neuroim2_exgridToIndex3DCpp(SEXP array_dimSEXP, SEXP iindSEXP, SEXP jindSEXP, SEXP kindSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type array_dim(array_dimSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type iind(iindSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type jind(jindSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type kind(kindSEXP);
    rcpp_result_gen = Rcpp::wrap(exgridToIndex3DCpp(array_dim, iind, jind, kind));
    return rcpp_result_gen;
END_RCPP
}
// exgridToIndex4DCpp
IntegerVector exgridToIndex4DCpp(IntegerVector array_dim, IntegerVector iind, IntegerVector jind, IntegerVector kind, IntegerVector mind);
RcppExport SEXP _neuroim2_exgridToIndex4DCpp(SEXP array_dimSEXP, SEXP iindSEXP, SEXP jindSEXP, SEXP kindSEXP, SEXP mindSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type array_dim(array_dimSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type iind(iindSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type jind(jindSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type kind(kindSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type mind(mindSEXP);
    rcpp_result_gen = Rcpp::wrap(exgridToIndex4DCpp(array_dim, iind, jind, kind, mind));
    return rcpp_result_gen;
END_RCPP
}
// gridToIndexCpp
IntegerVector gridToIndexCpp(IntegerVector array_dim, IntegerMatrix voxmat);
RcppExport SEXP _neuroim2_gridToIndexCpp(SEXP array_dimSEXP, SEXP voxmatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type array_dim(array_dimSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type voxmat(voxmatSEXP);
    rcpp_result_gen = Rcpp::wrap(gridToIndexCpp(array_dim, voxmat));
    return rcpp_result_gen;
END_RCPP
}
// gridToIndex3DCpp
IntegerVector gridToIndex3DCpp(IntegerVector array_dim, NumericMatrix voxmat);
RcppExport SEXP _neuroim2_gridToIndex3DCpp(SEXP array_dimSEXP, SEXP voxmatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type array_dim(array_dimSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type voxmat(voxmatSEXP);
    rcpp_result_gen = Rcpp::wrap(gridToIndex3DCpp(array_dim, voxmat));
    return rcpp_result_gen;
END_RCPP
}
// indexToGridCpp
NumericMatrix indexToGridCpp(IntegerVector idx, IntegerVector array_dim);
RcppExport SEXP _neuroim2_indexToGridCpp(SEXP idxSEXP, SEXP array_dimSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type idx(idxSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type array_dim(array_dimSEXP);
    rcpp_result_gen = Rcpp::wrap(indexToGridCpp(idx, array_dim));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_neuroim2_find_seqnum", (DL_FUNC) &_neuroim2_find_seqnum, 2},
    {"_neuroim2_grid_to_intvec", (DL_FUNC) &_neuroim2_grid_to_intvec, 2},
    {"_neuroim2_exgridToIndex3DCpp", (DL_FUNC) &_neuroim2_exgridToIndex3DCpp, 4},
    {"_neuroim2_exgridToIndex4DCpp", (DL_FUNC) &_neuroim2_exgridToIndex4DCpp, 5},
    {"_neuroim2_gridToIndexCpp", (DL_FUNC) &_neuroim2_gridToIndexCpp, 2},
    {"_neuroim2_gridToIndex3DCpp", (DL_FUNC) &_neuroim2_gridToIndex3DCpp, 2},
    {"_neuroim2_indexToGridCpp", (DL_FUNC) &_neuroim2_indexToGridCpp, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_neuroim2(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
