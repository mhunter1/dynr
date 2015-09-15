

dynr.run <- function(model,data) {
  tmp <- .Call("main_R", model, data, PACKAGE = "dynr")
  return(tmp)
}



