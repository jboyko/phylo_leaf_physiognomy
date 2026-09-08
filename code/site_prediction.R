# MAT is averaged in degrees C; log(MAP) is averaged on its fitted scale.
# Exponentiating the latter gives the geometric site mean used for fossils.
site_prediction <- function(predictions, target = c("mat", "log_map")) {
  target <- match.arg(target)
  predictions <- predictions[!is.na(predictions)]
  if (!length(predictions)) return(NA_real_)
  if (any(!is.finite(predictions))) stop("Non-finite species prediction")
  mean(predictions)
}
