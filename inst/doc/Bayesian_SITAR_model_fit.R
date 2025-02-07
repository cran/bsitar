## ----SETTINGS-knitr, include=FALSE------------------------------------------------------
stopifnot(require(knitr))
options(width = 90)
knitr::knit_hooks$set(purl = knitr::hook_purl)
opts_chunk$set(
  comment = NA,
  message = FALSE,
  warning = FALSE,
  # eval = if (isTRUE(exists("params"))) params$EVAL else FALSE,
  purl = FALSE,
  dev = "jpeg",
  dpi = 100,
  fig.cap = "",
  fig.asp = 0.8,
  fig.height=4,
  fig.width = 5,
  out.width = "60%",
  fig.align = "center"
)

