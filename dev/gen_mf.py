# Generates the .a/.r/.u.yaml files of psmodel and wtest. Run from the module root: python dev/gen_mf.py
import os

NL = chr(10)


def var(name, title, multi=False, perm=("numeric",), sug=None):
    s = f"    - name: {name}{NL}      title: {title}{NL}      type: {'Variables' if multi else 'Variable'}{NL}"
    if sug:
        s += f"      suggested:{NL}" + "".join(f"          - {x}{NL}" for x in sug)
    s += f"      permitted:{NL}" + "".join(f"          - {x}{NL}" for x in perm) + NL
    return s


def num(name, title, default, mn=None, mx=None, typ="Number"):
    s = f"    - name: {name}{NL}      title: {title}{NL}      type: {typ}{NL}      default: {default}{NL}"
    if mn is not None:
        s += f"      min: {mn}{NL}"
    if mx is not None:
        s += f"      max: {mx}{NL}"
    return s + NL


def lst(name, title, opts, default):
    s = f"    - name: {name}{NL}      title: {title}{NL}      type: List{NL}      options:{NL}"
    for n, t in opts:
        s += f"          - name: {n}{NL}            title: '{t}'{NL}"
    return s + f"      default: {default}{NL}{NL}"


def boo(name, title, default="false"):
    return f"    - name: {name}{NL}      title: {title}{NL}      type: Bool{NL}      default: {default}{NL}{NL}"


def head_a(name, title, sub, mtitle, msub):
    return (f"---{NL}name: {name}{NL}title: {title}{NL}jas: '1.2'{NL}version: '1.0.0'{NL}menuGroup: Matching{NL}menuSubgroup: {sub}{NL}"
            f"menuTitle: {mtitle}{NL}menuSubtitle: {msub}{NL}{NL}options:{NL}    - name: data{NL}      type: Data{NL}{NL}")


def tail_a():
    return boo("showHelp", "Method notes") + boo("showGuide", "Data and setup guide") + f"...{NL}"


def r_head(name, title):
    return (f"---{NL}name: {name}{NL}title: {title}{NL}jrs: '1.1'{NL}{NL}items:{NL}"
            f"    - name: guide{NL}      title: Data and setup guide{NL}      type: Html{NL}      visible: (showGuide){NL}{NL}"
            f"    - name: helpText{NL}      title: Method notes{NL}      type: Html{NL}      visible: (showHelp){NL}{NL}"
            f"    - name: notes{NL}      title: Notes{NL}      type: Html{NL}      visible: true{NL}{NL}")


def table(name, title, cols, clear, visible="true"):
    s = f"    - name: {name}{NL}      title: '{title}'{NL}      type: Table{NL}      visible: {visible}{NL}      clearWith:{NL}"
    for c in clear:
        s += f"          - {c}{NL}"
    s += f"      columns:{NL}"
    for n, t, ty in cols:
        s += f"          - name: {n}{NL}            title: '{t}'{NL}            type: {ty}{NL}"
        if n == "p":
            s += f"            format: zto,pvalue{NL}"
    return s + NL


def image(name, title, fun, clear, visible, w=640, h=400):
    s = (f"    - name: {name}{NL}      title: '{title}'{NL}      type: Image{NL}      width: {w}{NL}      height: {h}{NL}"
         f"      renderFun: {fun}{NL}      visible: {visible}{NL}      clearWith:{NL}")
    return s + "".join(f"          - {c}{NL}" for c in clear) + NL


def box(label, name, mx=None):
    s = (f"          - type: TargetLayoutBox{NL}            label: {label}{NL}            children:{NL}"
         f"                - type: VariablesListBox{NL}                  name: {name}{NL}                  isTarget: true{NL}"
         f"                  height: small{NL}")
    if mx:
        s += f"                  maxItemCount: {mx}{NL}"
    return s


def u_head(title, name, boxes):
    s = (f"title: {title}{NL}name: {name}{NL}jus: '3.0'{NL}stage: 0{NL}compilerMode: tame{NL}children:{NL}"
         f"    - type: VariableSupplier{NL}      persistentItems: false{NL}      stretchFactor: 1{NL}      children:{NL}")
    for lab, nm, mx in boxes:
        s += box(lab, nm, mx)
    return s


def collapse(label, children, collapsed=True):
    s = f"    - type: CollapseBox{NL}      label: {label}{NL}      collapsed: {'true' if collapsed else 'false'}{NL}      children:{NL}"
    for c in children:
        s += c
    return s


def combo(name):
    return f"          - type: ComboBox{NL}            name: {name}{NL}"


def check(name, label=None):
    s = f"          - type: CheckBox{NL}            name: {name}{NL}"
    if label:
        s += f"            label: {label}{NL}"
    return s


def text(name, label, enable=None):
    s = f"          - type: TextBox{NL}            name: {name}{NL}            label: {label}{NL}            format: number{NL}"
    if enable:
        s += f"            enable: {enable}{NL}"
    return s


def write(path, s):
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write(s)


OUT = collapse("Output", [check("showHelp", "Method notes"), check("showGuide", "Data and setup guide")])
TRIM = (lst("trim", "Trim observations by propensity score", [("none", "No trimming"), ("minmax", "Outside the common range of the two groups"),
                                                                   ("crump", "Outside 0.1 to 0.9"), ("custom", "Outside the bounds below")], "none") +
        num("trimLo", "Lower bound", 0.05, 0, 0.5) + num("trimHi", "Upper bound", 0.95, 0.5, 1))
LINK = lst("link", "Propensity score model", [("logit", "Logit"), ("probit", "Probit")], "logit")
LEVEL = num("level", "Confidence level (%)", 95, 50, 99.9)

os.makedirs("jamovi", exist_ok=True)

# ---------------------------------------------------------------- psmodel
write("jamovi/psmodel.a.yaml",
      head_a("psmodel", "Propensity Score and Overlap", "Preparation", "Propensity Score and Overlap", "Model, overlap and trimming") +
      var("treat", "Treatment (0/1)", perm=("numeric", "factor")) +
      var("covs", "Covariates", multi=True, perm=("numeric", "factor"), sug=("continuous", "nominal")) +
      var("weights", "Sampling weights") + LINK + TRIM + LEVEL + boo("plot", "Propensity score distribution", "true") + tail_a())
write("jamovi/psmodel.r.yaml",
      r_head("psmodel", "Propensity Score and Overlap") +
      table("coef", "Propensity score model", [("term", "Variable", "text"), ("estimate", "Coefficient", "number"), ("se", "Std. error", "number"),
                                              ("z", "z", "number"), ("p", "p", "number"), ("lo", "CI lower", "number"), ("hi", "CI upper", "number")],
            ["treat", "covs", "weights", "link", "level"]) +
      table("overlap", "Propensity score by group", [("group", "Group", "text"), ("n", "N", "integer"), ("min", "Min", "number"), ("q05", "5th percentile", "number"),
                                                    ("median", "Median", "number"), ("q95", "95th percentile", "number"), ("max", "Max", "number")],
            ["treat", "covs", "weights", "link"]) +
      table("summary", "Model summary and common support", [("item", "", "text"), ("value", "", "text")],
            ["treat", "covs", "weights", "link", "trim", "trimLo", "trimHi"]) +
      image("psPlot", "Propensity score distribution by treatment group", ".plotPs", ["treat", "covs", "weights", "link", "trim", "trimLo", "trimHi"], "(plot)") +
      f"...{NL}")
write("jamovi/psmodel.u.yaml",
      u_head("Propensity Score and Overlap", "psmodel", [("Treatment (0/1)", "treat", 1), ("Covariates", "covs", None), ("Sampling weights", "weights", 1)]) +
      collapse("Model and trimming", [combo("link"), combo("trim"), text("trimLo", "Lower bound", "(trim:custom)"), text("trimHi", "Upper bound", "(trim:custom)"),
                                      text("level", "Confidence level (%)"), check("plot")]) + OUT)

# ---------------------------------------------------------------- wtest
write("jamovi/wtest.a.yaml",
      head_a("wtest", "Weighting and Doubly Robust Estimators", "Estimation", "Weighting and Doubly Robust Estimators", "RA, IPW, AIPW, IPWRA") +
      var("outcome", "Outcome", sug=("continuous",)) + var("treat", "Treatment (0/1)", perm=("numeric", "factor")) +
      var("covs", "Covariates", multi=True, perm=("numeric", "factor"), sug=("continuous", "nominal")) +
      var("extra", "Additional outcome-model covariates", multi=True, perm=("numeric", "factor"), sug=("continuous", "nominal")) +
      var("weights", "Sampling weights") + var("cluster", "Cluster variable", perm=("factor", "id", "numeric")) +
      lst("estimand", "Estimand", [("all", "All three: ATE, ATT and ATC"), ("ate", "Average treatment effect (ATE)"), ("att", "Effect on the treated (ATT)"), ("atc", "Effect on the controls (ATC)")], "all") +
      boo("doRa", "Regression adjustment (RA)", "true") + boo("doIpw", "Inverse probability weighting (IPW)", "true") +
      boo("doAipw", "Augmented IPW (AIPW)", "true") + boo("doIpwra", "IPW regression adjustment (IPWRA)", "true") +
      LINK + TRIM + LEVEL + boo("showPom", "Potential-outcome means", "false") + boo("showWeights", "Weight diagnostics", "true") + tail_a())
write("jamovi/wtest.r.yaml",
      r_head("wtest", "Weighting and Doubly Robust Estimators") +
      table("effects", "Treatment effect estimates", [("estimand", "Estimand", "text"), ("method", "Estimator", "text"), ("estimate", "Effect", "number"), ("se", "Std. error", "number"),
                                                     ("z", "z", "number"), ("p", "p", "number"), ("lo", "CI lower", "number"), ("hi", "CI upper", "number")],
            ["outcome", "treat", "covs", "extra", "weights", "cluster", "estimand", "doRa", "doIpw", "doAipw", "doIpwra", "link", "trim", "trimLo", "trimHi", "level"]) +
      table("pom", "Potential-outcome means", [("estimand", "Estimand", "text"), ("method", "Estimator", "text"), ("group", "Mean of", "text"), ("estimate", "Estimate", "number"), ("se", "Std. error", "number")],
            ["outcome", "treat", "covs", "extra", "weights", "cluster", "estimand", "doRa", "doIpw", "doAipw", "doIpwra", "link", "trim", "trimLo", "trimHi"], "(showPom)") +
      table("wdiag", "Inverse-probability weights", [("estimand", "Estimand", "text"), ("group", "Group", "text"), ("n", "N", "integer"), ("min", "Min", "number"), ("mean", "Mean", "number"),
                                                     ("max", "Max", "number"), ("ess", "Effective N", "number"), ("share", "Largest weight, share of the total", "number")],
            ["treat", "covs", "weights", "estimand", "link", "trim", "trimLo", "trimHi"], "(showWeights)") +
      f"...{NL}")
write("jamovi/wtest.u.yaml",
      u_head("Weighting and Doubly Robust Estimators", "wtest",
             [("Outcome", "outcome", 1), ("Treatment (0/1)", "treat", 1), ("Covariates (propensity score)", "covs", None),
              ("Extra outcome-model covariates", "extra", None), ("Sampling weights", "weights", 1), ("Cluster variable", "cluster", 1)]) +
      collapse("Estimand and estimators", [combo("estimand"), check("doRa"), check("doIpw"), check("doAipw"), check("doIpwra"), combo("link"), text("level", "Confidence level (%)")]) +
      collapse("Trimming", [combo("trim"), text("trimLo", "Lower bound", "(trim:custom)"), text("trimHi", "Upper bound", "(trim:custom)")]) +
      collapse("Tables", [check("showPom"), check("showWeights")]) + OUT)
print("written")
