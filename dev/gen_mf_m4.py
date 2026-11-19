# Generates the .a/.r/.u.yaml files of matchcmp and sensobs (M4). Run from the module root: python dev/gen_mf_m4.py
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
from gen_mf_m3 import *  # noqa  (re-writes the earlier analyses, which is harmless)

# ---------------------------------------------------------------- matchcmp
MC_CLEAR = ["outcome", "treat", "covs", "estimand", "link", "k", "biasadj", "moments", "doNaive", "doRa", "doIpw", "doAipw", "doIpwra", "doNn", "doEbal"]
write("jamovi/matchcmp.a.yaml",
      head_a("matchcmp", "Compare Matching and Weighting Estimators", "Estimation", "Compare Estimators", "Side by side, same data and estimand") +
      var("outcome", "Outcome", sug=("continuous",)) + var("treat", "Treatment (0/1)", perm=("numeric", "factor")) +
      var("covs", "Covariates", multi=True, perm=("numeric", "factor"), sug=("continuous", "nominal")) +
      ESTIMAND_ALL + LINK +
      boo("doNaive", "Unadjusted difference in means", "true") + boo("doRa", "Regression adjustment (RA)", "true") + boo("doIpw", "Inverse probability weighting (IPW)", "true") +
      boo("doAipw", "Augmented IPW (AIPW)", "true") + boo("doIpwra", "IPW regression adjustment (IPWRA)", "true") +
      boo("doNn", "Nearest-neighbour matching", "true") + boo("doEbal", "Entropy balancing", "true") +
      num("k", "Number of neighbours for matching", 1, 1, 50, "Integer") + boo("biasadj", "Bias adjustment for matching", "true") +
      lst("moments", "Moments for entropy balancing", [("m1", "Means"), ("m2", "Means and variances"), ("m3", "Means, variances and skewness")], "m1") +
      LEVEL + boo("plot", "Forest plot", "true") + tail_a())
write("jamovi/matchcmp.r.yaml",
      r_head("matchcmp", "Compare Matching and Weighting Estimators") +
      table2("est", "Estimates by method", [("estimand", "Estimand", "text"), ("method", "Estimator", "text"), ("estimate", "Effect", "number"), ("se", "Std. error", "number"), ("p", "p", "number"),
                                            ("lo", "CI lower", "number"), ("hi", "CI upper", "number"),
                                            ("smd", "Largest |std. diff.|", "number", "Balance after adjustment"), ("above", "Variables above 0.1", "integer", "Balance after adjustment")],
             MC_CLEAR + ["level"]) +
      image("forest", "Effect estimates with confidence intervals", ".plotForest", MC_CLEAR + ["level"], "(plot)", 640, 560) +
      f"...{NL}")
write("jamovi/matchcmp.u.yaml",
      u_head("Compare Matching and Weighting Estimators", "matchcmp",
             [("Outcome", "outcome", 1), ("Treatment (0/1)", "treat", 1), ("Covariates", "covs", None)]) +
      collapse("Estimand and estimators", [combo("estimand"), combo("link"), check("doNaive"), check("doRa"), check("doIpw"), check("doAipw"), check("doIpwra"), check("doNn"),
                                           check("doEbal"), text("level", "Confidence level (%)")]) +
      collapse("Matching and balancing options", [text("k", "Number of neighbours for matching"), check("biasadj"), combo("moments")], True) +
      collapse("Plot", [check("plot")]) + OUT)

# ---------------------------------------------------------------- sensobs
SO_CLEAR = ["mode", "outcome", "treat", "covs", "estimand", "link", "mType", "mEst", "mLo", "mHi", "k", "level"]
write("jamovi/sensobs.a.yaml",
      head_a("sensobs", "Sensitivity to Unobserved Confounding", "Sensitivity", "Sensitivity to Unobserved Confounding", "E-value and Rosenbaum bounds") +
      lst("mode", "Effect to assess", [("data", "Estimate it from the data (AIPW)"), ("manual", "Enter an estimate and its confidence limits")], "data") +
      var("outcome", "Outcome", sug=("continuous",)) + var("treat", "Treatment (0/1)", perm=("numeric", "factor")) +
      var("covs", "Covariates", multi=True, perm=("numeric", "factor"), sug=("continuous", "nominal")) +
      ESTIMAND + LINK +
      lst("mType", "Effect measure", [("smd", "Standardised mean difference"), ("rr", "Risk ratio"), ("or_rare", "Odds ratio, rare outcome"),
                                      ("or_common", "Odds ratio, common outcome"), ("hr_rare", "Hazard ratio, rare outcome"), ("hr_common", "Hazard ratio, common outcome")], "rr") +
      num("mEst", "Estimate", 1.5) + num("mLo", "Lower confidence limit", 1.2) + num("mHi", "Upper confidence limit", 1.9) +
      boo("doRosenbaum", "Rosenbaum bounds (matched pairs)", "true") +
      num("k", "Matches per treated unit (Rosenbaum)", 1, 1, 10, "Integer") +
      num("gammaMax", "Largest Gamma (Rosenbaum)", 3, 1.1, 20) + num("gammaStep", "Gamma step (Rosenbaum)", 0.25, 0.01, 5) +
      LEVEL + boo("plot", "Plot of the p-value bounds", "true") + tail_a())
write("jamovi/sensobs.r.yaml",
      r_head("sensobs", "Sensitivity to Unobserved Confounding") +
      table("effect", "Estimate assessed", [("method", "Estimand", "text"), ("estimate", "Effect", "number"), ("se", "Std. error", "number"), ("lo", "CI lower", "number"),
                                            ("hi", "CI upper", "number"), ("d", "Std. mean difference", "number"), ("dlo", "SMD lower", "number"), ("dhi", "SMD upper", "number")],
            SO_CLEAR, "(mode:data)") +
      table("evalue", "E-value", [("item", "", "text"), ("value", "Value", "number"), ("rr", "Risk-ratio equivalent", "number"), ("e", "E-value", "number")], SO_CLEAR) +
      table("rsum", "Rosenbaum bounds: matched pairs", [("item", "", "text"), ("value", "", "text")],
            ["mode", "outcome", "treat", "covs", "k", "doRosenbaum"], "(mode:data && doRosenbaum)") +
      table("rosen", "Rosenbaum bounds: one-sided p-values", [("gamma", "Gamma", "number"), ("wlo", "Signed rank, lower bound", "number"), ("whi", "Signed rank, upper bound", "number"),
                                                              ("slo", "Sign test, lower bound", "number"), ("shi", "Sign test, upper bound", "number")],
            ["mode", "outcome", "treat", "covs", "k", "doRosenbaum", "gammaMax", "gammaStep"], "(mode:data && doRosenbaum)") +
      image("rosPlot", "Upper bound of the p-value as hidden bias grows", ".plotRos", ["mode", "outcome", "treat", "covs", "k", "doRosenbaum", "gammaMax", "gammaStep"],
            "(mode:data && doRosenbaum && plot)", 560, 360) +
      f"...{NL}")
write("jamovi/sensobs.u.yaml",
      u_head("Sensitivity to Unobserved Confounding", "sensobs",
             [("Outcome", "outcome", 1), ("Treatment (0/1)", "treat", 1), ("Covariates", "covs", None)]) +
      collapse("Effect to assess", [combo("mode"), combo("estimand"), combo("link"), text("level", "Confidence level (%)")], False) +
      collapse("Manual input", [combo("mType"), text("mEst", "Estimate", "(mode:manual)"), text("mLo", "Lower confidence limit", "(mode:manual)"),
                                text("mHi", "Upper confidence limit", "(mode:manual)")], True) +
      collapse("Rosenbaum bounds", [check("doRosenbaum"), text("k", "Matches per treated unit"), text("gammaMax", "Largest Gamma"), text("gammaStep", "Step"), check("plot")], True) +
      OUT)

# ---------------------------------------------------------------- references: added to each .a.yaml after it is written
REFS = {"psmodel": ["rosenbaum1983", "crump2009"], "wtest": ["rosenbaum1983", "wooldridge2007", "robins1994", "imbensrubin2015"],
        "balance": ["austin2009", "hainmueller2012"], "ebal": ["hainmueller2012"], "nnmatch": ["abadieimbens2006", "abadieimbens2011"],
        "matchcmp": ["imbensrubin2015"], "sensobs": ["vanderweeleding2017", "rosenbaum2002"]}
for name, ids in REFS.items():
    path = f"jamovi/{name}.r.yaml"
    txt = open(path, encoding="utf-8").read()
    if NL + "refs:" not in txt:
        txt = txt.rstrip()
        assert txt.endswith("...")
        txt = txt[:-3] + "refs:" + NL + "".join(f"    - {i}{NL}" for i in ids) + NL + "..." + NL
        write(path, txt)
print("written m4")
