library("devtools")
pkg_path = '/home/pals2/Work/ggseqlogo'

devtools::build(pkg_path)
devtools::install(pkg_path)
