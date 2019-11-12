
.onAttach <- function(lib, pkg)  {
  packageStartupMessage("version: ", utils::packageDescription("landmap", field="Version"), appendLF = TRUE)
}
