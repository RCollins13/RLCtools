# RLCtools
## An R library of helper functions  

Copyright (c) 2023-Present, [Ryan L. Collins](mailto:Ryan_Collins@dfci.harvard.edu) and the Dana-Farber Cancer Institute.  
Distributed under terms of the [GNU GPL v2.0 License](/LICENSE) (see `LICENSE`).  

---  

`RLCtools` is an R library of helper functions used across projects.  

It can be installed from source like any other R package. For example:  
```
version <- 0.1
install.packages(paste("RLCtools_", version, ".tar.gz", sep=""), repos=NULL, type="source")
```

Once installed, you load it like any other R package:
```
library(RLCtools)
```

Most functions have help text, which can be accessed using `?` like any other R function.  

There is a companion Docker image with minimal libraries beyond those required for `RLCtools`. See the `docker` subdirectory for more information.  
