# Packages the pipeline depends on at runtime but never names, so renv cannot find
# them by reading the scripts. This file exists only to declare them; nothing sources it.
#
# ggplot2::ggsave() picks the PNG device by availability: ragg::agg_png() when ragg is
# installed, grDevices::png() when it is not. The two shape text differently, so a
# library without ragg silently renders every committed figure to different bytes --
# no error, no warning. ragg is load-bearing for the figures even though no script
# calls it. textshaping comes with it.
library(ragg)
