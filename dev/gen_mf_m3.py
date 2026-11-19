# Generates the .a/.r/.u.yaml files of nnmatch (M3). Run from the module root: python dev/gen_mf_m3.py
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
from gen_mf_m2 import *  # noqa  (re-writes the earlier analyses, which is harmless)

NN_CLEAR = ["outcome", "treat", "covs", "adjvars", "exact", "estimand", "k", "metric", "biasadj", "caliper", "varNn"]
write("jamovi/nnmatch.a.yaml",
      head_a("nnmatch", "Matching Estimators", "Estimation", "Matching Estimators", "Nearest-neighbour matching") +
      var("outcome", "Outcome", sug=("continuous",)) + var("treat", "Treatment (0/1)", perm=("numeric", "factor")) +
      var("covs", "Matching covariates", multi=True, perm=("numeric", "factor"), sug=("continuous", "nominal")) +
      var("adjvars", "Bias-adjustment covariates", multi=True, perm=("numeric", "factor"), sug=("continuous", "nominal")) +
      var("exact", "Exact-match variables", multi=True, perm=("numeric", "factor"), sug=("nominal", "ordinal")) +
      ESTIMAND_ALL +
      num("k", "Number of neighbours", 1, 1, 50, "Integer") +
      lst("metric", "Distance", [("maha", "Mahalanobis"), ("ivar", "Inverse variance (diagonal)"), ("eucl", "Euclidean")], "maha") +
      boo("biasadj", "Bias adjustment (Abadie and Imbens 2011)", "false") +
      num("caliper", "Caliper (distance units; 0 = none)", 0, 0) +
      num("varNn", "Neighbours for the variance estimate", 2, 1, 20, "Integer") +
      LEVEL + boo("showMatches", "Matches per unit", "true") + boo("showBalance", "Balance after matching", "true") + tail_a())
write("jamovi/nnmatch.r.yaml",
      r_head("nnmatch", "Matching Estimators") +
      table("effects", "Treatment effect by nearest-neighbour matching", [("method", "Estimand", "text"), ("estimate", "Effect", "number"), ("se", "Std. error", "number"),
                                                                           ("z", "z", "number"), ("p", "p", "number"), ("lo", "CI lower", "number"), ("hi", "CI upper", "number")],
            NN_CLEAR + ["level"]) +
      table("summary", "Matching summary", [("estimand", "Estimand", "text"), ("item", "", "text"), ("value", "", "text")], NN_CLEAR) +
      table("uses", "Times a unit is used as a match", [("estimand", "Estimand", "text"), ("group", "Group", "text"), ("n0", "0", "integer"), ("n1", "1", "integer"), ("n2", "2", "integer"),
                                                       ("n3", "3 to 4", "integer"), ("n5", "5 or more", "integer"), ("max", "Maximum", "number")],
            NN_CLEAR, "(showMatches)") +
      table2("bal", "Covariate balance", [("estimand", "Estimand", "text"), ("variable", "Variable", "text"),
                                          ("m1", "Treated mean", "number", "Unadjusted"), ("m0", "Control mean", "number", "Unadjusted"), ("smd", "Std. diff.", "number", "Unadjusted"),
                                          ("am1", "Treated mean", "number", "Matched"), ("am0", "Control mean", "number", "Matched"),
                                          ("asmd", "Std. diff.", "number", "Matched"), ("avr", "Var. ratio", "number", "Matched")],
             NN_CLEAR, "(showBalance)") +
      f"...{NL}")
write("jamovi/nnmatch.u.yaml",
      u_head("Matching Estimators", "nnmatch",
             [("Outcome", "outcome", 1), ("Treatment (0/1)", "treat", 1), ("Matching covariates", "covs", None),
              ("Bias-adjustment covariates (default matching covariates)", "adjvars", None), ("Exact-match variables", "exact", None)]) +
      collapse("Estimand and matching", [combo("estimand"), text("k", "Number of neighbours"), combo("metric"), text("caliper", "Caliper (distance units; 0 = none)"),
                                         check("biasadj"), text("varNn", "Neighbours for the variance estimate"), text("level", "Confidence level (%)")]) +
      collapse("Tables", [check("showMatches"), check("showBalance")]) + OUT)
print("written m3")
