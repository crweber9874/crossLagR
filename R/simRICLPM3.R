#' @title simRICLPM
#' @description Simulate data from an observed, random intercept cross-lagged regression (RI-CLR) model
#' In particular, this creates a synthetic dataset of a cross-lagged model with a specified number of waves and structural parameters.
#'
#' @param waves The number of waves (time points) in the model.
#' @param stability.p The stability parameter for the x variable (autoregressive effect).
#' @param stability.q The stability parameter for the y variable (autoregressive effect).
#' @param stability.r The stability parameter for the z variable (autoregressive effect).
#' @param cov.pq The covariance between x and y within the same time point.
#' @param cov.pr The covariance between x and z within the same time point.
#' @param cov.qr The covariance between y and z within the same time point.
#' @param cross.q The cross-lagged effect of x on y at the next time point.
#' @param cross.p The cross-lagged effect of y on x at the next time point.
#' @param cross.r The cross-lagged effect of z on x and y at the next time point.
#' @param ... Additional arguments to pass to the `lavaan::simulateData` function.
#'
#' @return A list containing two elements:
#'    * `model`: The Lavaan model syntax used for data simulation.
#'    * `data`:  The simulated data in a data frame format.
#'
#' @export
#'

simRICLPM3 <- function(waves = 10,
                       stability.p = 0.2,
                       stability.q = 0.2,
                       stability.r = 0.2,
                       cross.p = 0.1,
                       cross.q = 0.1,
                       cross.r = 0.1,
                       variance.p = 1,
                       variance.q = 1,
                       variance.r = 1,
                       cov.pq = 0.1,
                       cov.pr = 0.1,
                       cov.qr = 0.1,
                       variance.between.x = 1, # random intercepts, x
                       variance.between.y = 1, # random intercepts, y
                       variance.between.z = 1, # random intercepts, z
                       cov.between.xy = 0.5, # covariance of intercept terms x and y
                       cov.between.xz = 0.5, # covariance of intercept terms x and z
                       cov.between.yz = 0.5, # covariance of intercept terms y and z
                       ...) {
        model_string <- ""
        model_string <- paste0(model_string, "\n BX =~  1* x1")
        for (w in 2:waves) {
                model_string <- paste0(model_string, " + 1 *x", w, "")
        }
        model_string <- paste0(model_string, "\n BY =~   1* y1")
        for (w in 2:waves) {
                model_string <- paste0(model_string, " + 1 * y", w)
        }
        model_string <- paste0(model_string, "\n BZ =~   1* z1")
        for (w in 2:waves) {
                model_string <- paste0(model_string, " + 1 * z", w)
        }
        model_string <- paste0(model_string, "\n")

        for (w in 1:waves) {
                model_string <- paste0(model_string, "x", w, "~ 1", "\n")
        }
        for (w in 1:waves) {
                model_string <- paste0(model_string, "y", w, "~ 1", "\n")
        }
        for (w in 1:waves) {
                model_string <- paste0(model_string, "z", w, "~ 1", "\n")
        }
        for (w in 1:waves) {
                model_string <- paste0(
                        model_string, "\np", w, " =~ 1*x", w,
                        "\nq", w, " =~ 1*y", w,
                        "\nr", w, " =~ 1*z", w
                )
        }

        # Stability
        for (w in 2:waves) {
                model_string <- paste0(
                        model_string, "\n p", w, " ~ ", stability.p, " * p", w - 1, " + ", cross.q, " * q", w - 1, " + ", cross.r, " * r", w - 1,
                        "\n q", w, " ~  ", stability.q, " * q", w - 1, " + ", cross.p, " * p", w - 1, " + ", cross.r, " * r", w - 1,
                        "\n r", w, " ~  ", stability.r, " * r", w - 1, " + ", cross.p, " * p", w - 1, " + ", cross.q, " * q", w - 1
                )
        }

        for (w in 1:waves) {
                model_string <- paste0(
                        model_string, "\n p", w, " ~~ ", variance.p, " * p", w,
                        "\n q", w, " ~~ ", variance.q, " * q", w,
                        "\n r", w, " ~~ ", variance.r, " * r", w,
                        "\n p", w, " ~~ ", cov.pq, " * q", w,
                        "\n p", w, " ~~ ", cov.pr, " * r", w,
                        "\n q", w, " ~~ ", cov.qr, " * r", w
                )
        }

        for (w in 1:waves) {
                model_string <- paste0(
                        model_string, "\nx", w, "~~",
                        "0*", "x", w
                )
        }
        for (w in 1:waves) {
                model_string <- paste0(
                        model_string, "\ny", w, "~~",
                        "0*", "y", w
                )
        }
        for (w in 1:waves) {
                model_string <- paste0(
                        model_string, "\nz", w, "~~",
                        "0*", "z", w
                )
        }

        model_string <- paste0(
                model_string,
                "\n BX ~~", variance.between.x, "* BX",
                "\n BY ~~", variance.between.y, "* BY",
                "\n BZ ~~", variance.between.z, "* BZ",
                "\n BX ~~", cov.between.xy, "* BY",
                "\n BX ~~", cov.between.xz, "* BZ",
                "\n BY ~~", cov.between.yz, "* BZ"
        )

        dat <- lavaan::simulateData(
                model = model_string,
                int.ov.free = FALSE
        )

        return(list(model = model_string, data = dat))
}
