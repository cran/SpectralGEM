.First.lib <- function(lib,pkg)
{
   library.dynam("SpectralGEM",pkg,lib)
   cat("SpectralGEM 0.1-1 loaded\n")
}

