# Generates the .a/.r/.u.yaml files of balance and ebal (M2). Run from the module root: python dev/gen_mf_m2.py
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
from gen_mf import *  # noqa  (re-writes psmodel and wtest, which is harmless)

ESTIMAND_ALL = lst("estimand", "Estimand", [("all", "All three: ATE, ATT and ATC"), ("ate", "Average treatment effect (ATE)"), ("att", "Effect on the treated (ATT)"), ("atc", "Effect on the controls (ATC)")], "all")
ESTIMAND = lst("estimand", "Estimand", [("ate", "Average treatment effect (ATE)"), ("att", "Effect on the treated (ATT)"), ("atc", "Effect on the controls (ATC)")], "ate")


def table2(name, title, cols, clear, visible="true"):
    """Like table() but a column may carry a superTitle (4th element) and a visibility rule (5th)."""
    s = f"    - name: {name}{NL}      title: '{title}'{NL}      type: Table{NL}      visible: {visible}{NL}      clearWith:{NL}"
    for c in clear:
        s += f"          - {c}{NL}"
    s += f"      columns:{NL}"
    for c in cols:
        n, t, ty = c[:3]
        s += f"          - name: {n}{NL}            title: '{t}'{NL}            type: {ty}{NL}"
        if len(c) > 3 and c[3]:
            s += f"            superTitle: {c[3]}{NL}"
        if len(c) > 4 and c[4]:
            s += f"            visible: {c[4]}{NL}"
        if n == "p":
            s += f"            format: zto,pvalue{NL}"
    return s + NL


# ---------------------------------------------------------------- balance
write("jamovi/balance.a.yaml",
      head_a("balance", "Covariate Balance", "Preparation", "Covariate Balance", "Standardised differences, before and after weighting") +
      var("treat", "Treatment (0/1)", perm=("numeric", "factor")) +
      var("covs", "Covariates", multi=True, perm=("numeric", "factor"), sug=("continuous", "nominal")) +
      var("weights", "Sampling weights") + var("adjW", "Weights variable") +
      lst("adjust", "Adjustment", [("none", "None (unadjusted only)"), ("ipw", "Inverse probability weights"), ("ebal", "Entropy balancing"),
                                   ("var", "Weights from a variable")], "ipw") +
      ESTIMAND + LINK + lst("moments", "Moments to balance (entropy balancing)", [("m1", "Means"), ("m2", "Means and variances"), ("m3", "Means, variances and skewness")], "m1") +
      lst("denom", "Standardise by", [("auto", "Follows the estimand"), ("pooled", "Pooled standard deviation"),
                                      ("treated", "Treated group"), ("control", "Control group")], "auto") +
      boo("stdBinary", "Standardise binary variables", "true") +
      num("threshold", "Balance threshold (absolute standardised difference)", 0.1, 0.001, 5) +
      boo("showMeans", "Group means", "true") + boo("showVr", "Variance ratio", "true") + boo("showKs", "Kolmogorov-Smirnov statistic", "true") +
      boo("plot", "Love plot", "true") + tail_a())
SMD_CLEAR = ["treat", "covs", "weights", "adjW", "adjust", "estimand", "link", "moments", "denom", "stdBinary"]
write("jamovi/balance.r.yaml",
      r_head("balance", "Covariate Balance") +
      table2("summary", "Balance summary", [("item", "", "text"), ("before", "Unadjusted", "text"), ("after", "Adjusted", "text", None, "(adjust:ipw || adjust:ebal || adjust:var)")],
             SMD_CLEAR + ["threshold"]) +
      table2("bal", "Covariate balance", [("variable", "Variable", "text"),
                                          ("m1", "Treated mean", "number", "Unadjusted", "(showMeans)"), ("m0", "Control mean", "number", "Unadjusted", "(showMeans)"),
                                          ("smd", "Std. diff.", "number", "Unadjusted"), ("vr", "Var. ratio", "number", "Unadjusted", "(showVr)"),
                                          ("ks", "KS", "number", "Unadjusted", "(showKs)"),
                                          ("am1", "Treated mean", "number", "Adjusted", "(showMeans && (adjust:ipw || adjust:ebal || adjust:var))"),
                                          ("am0", "Control mean", "number", "Adjusted", "(showMeans && (adjust:ipw || adjust:ebal || adjust:var))"),
                                          ("asmd", "Std. diff.", "number", "Adjusted", "(adjust:ipw || adjust:ebal || adjust:var)"),
                                          ("avr", "Var. ratio", "number", "Adjusted", "(showVr && (adjust:ipw || adjust:ebal || adjust:var))"),
                                          ("aks", "KS", "number", "Adjusted", "(showKs && (adjust:ipw || adjust:ebal || adjust:var))")],
             SMD_CLEAR) +
      image("lovePlot", "Absolute standardised differences", ".plotLove", SMD_CLEAR + ["threshold"], "(plot)", 640, 420) +
      f"...{NL}")
write("jamovi/balance.u.yaml",
      u_head("Covariate Balance", "balance",
             [("Treatment (0/1)", "treat", 1), ("Covariates", "covs", None), ("Sampling weights", "weights", 1), ("Weights variable (for 'Weights from a variable')", "adjW", 1)]) +
      collapse("Adjustment", [combo("adjust"), combo("estimand"), combo("link"), combo("moments")]) +
      collapse("Statistics", [combo("denom"), check("stdBinary"), text("threshold", "Balance threshold"), check("showMeans"), check("showVr"), check("showKs"), check("plot")]) +
      OUT)

# ---------------------------------------------------------------- ebal
EB_CLEAR = ["outcome", "treat", "covs", "weights", "cluster", "estimand", "moments"]
write("jamovi/ebal.a.yaml",
      head_a("ebal", "Entropy Balancing", "Estimation", "Entropy Balancing", "Weights that match the covariate moments exactly") +
      var("outcome", "Outcome", sug=("continuous",)) + var("treat", "Treatment (0/1)", perm=("numeric", "factor")) +
      var("covs", "Covariates", multi=True, perm=("numeric", "factor"), sug=("continuous", "nominal")) +
      var("weights", "Sampling weights") + var("cluster", "Cluster variable", perm=("factor", "id", "numeric")) +
      ESTIMAND_ALL + lst("moments", "Moments to balance", [("m1", "Means"), ("m2", "Means and variances"), ("m3", "Means, variances and skewness")], "m1") +
      LEVEL + boo("showPom", "Potential-outcome means", "false") + boo("showWeights", "Weight diagnostics", "true") +
      boo("showBalance", "Balance table", "true") + tail_a())
write("jamovi/ebal.r.yaml",
      r_head("ebal", "Entropy Balancing") +
      table("effects", "Treatment effect after entropy balancing", [("method", "Estimand", "text"), ("estimate", "Effect", "number"), ("se", "Std. error", "number"),
                                                                   ("z", "z", "number"), ("p", "p", "number"), ("lo", "CI lower", "number"), ("hi", "CI upper", "number")],
            EB_CLEAR + ["level"]) +
      table("pom", "Potential-outcome means", [("estimand", "Estimand", "text"), ("group", "Mean of", "text"), ("estimate", "Estimate", "number"), ("se", "Std. error", "number")], EB_CLEAR, "(showPom)") +
      table("wdiag", "Entropy balancing weights", [("estimand", "Estimand", "text"), ("group", "Group", "text"), ("n", "N", "integer"), ("min", "Min", "number"), ("mean", "Mean", "number"),
                                                   ("max", "Max", "number"), ("ess", "Effective N", "number"), ("share", "Largest weight, share of the total", "number")],
            EB_CLEAR[1:], "(showWeights)") +
      table2("bal", "Covariate balance", [("estimand", "Estimand", "text"), ("variable", "Variable", "text"),
                                          ("m1", "Treated mean", "number", "Unadjusted"), ("m0", "Control mean", "number", "Unadjusted"), ("smd", "Std. diff.", "number", "Unadjusted"),
                                          ("am1", "Treated mean", "number", "After balancing"), ("am0", "Control mean", "number", "After balancing"),
                                          ("asmd", "Std. diff.", "number", "After balancing"),
                                          ("avr", "Var. ratio", "number", "After balancing")],
             EB_CLEAR[1:], "(showBalance)") +
      f"...{NL}")
write("jamovi/ebal.u.yaml",
      u_head("Entropy Balancing", "ebal",
             [("Outcome", "outcome", 1), ("Treatment (0/1)", "treat", 1), ("Covariates", "covs", None), ("Sampling weights", "weights", 1), ("Cluster variable", "cluster", 1)]) +
      collapse("Estimand and moments", [combo("estimand"), combo("moments"), text("level", "Confidence level (%)")]) +
      collapse("Tables", [check("showPom"), check("showWeights"), check("showBalance")]) + OUT)
print("written m2")
