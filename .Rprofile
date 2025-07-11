source("renv/activate.R")
renv::autoload()
if (dir.exists(".git")) {
  try({
    message("Running git pull...")
    result <- if (.Platform$OS.type == "windows") {
      shell("git pull", intern = TRUE)
    } else {
      system("git pull", intern = TRUE)
    }
    cat(result, sep = "\n")
  }, silent = TRUE)
}
