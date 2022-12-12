## Great-circle_distance
.dot <- function(x, y) {
  if (is.vector(x) && is.vector(y)) {
    sum(x * y)
  } else if (is.vector(x) && is.matrix(y)) {
    y %*% as.matrix(x)
  } else if (is.vector(y) && is.matrix(x)) {
    x %*% as.matrix(y)
  } else if (is.matrix(x) && is.matrix(y)) {
    apply(x * y, 1, sum)
  }
}

.norm <- function(x) {
  if (is.vector(x)) {
    sqrt(sum(x^2))
  } else if (is.matrix(x)) {
    apply(x, 1, function(x) sqrt(sum(x^2)))
  }
}

is_coord <- function(x) {
  if (is.vector(x)) {
    length(x) == 2
  } else if (is.matrix(x)) {
    ncol(x) == 2
  } else {
    FALSE
  }
}

sel_x <- function(x) {
  if (is.vector(x)) {
    x[1]
  } else if (is.matrix(x)) {
    x[, 1]
  }
}

sel_y <- function(x) {
  if (is.vector(x)) {
    x[2]
  } else if (is.matrix(x)) {
    x[, 2]
  }
}

#' Rotation angle of two speed vector
#'
#' @param v1,v2 speed vector, in deg, (lon, lat).
#' - `v1` is the previous time.
#' - `v2` is the current time;
#' @param clockwise one of `c(-1, 1)`:
#' -`clockwise = -1`: clockwise corresponds to negative angle
#' -`clockwise =  1`: clockwise corresponds to positive angle (azimuth angle in this style)
#' 
#' @return
#' **v2 -> v1**:
#' 
#' * rotation angle (in deg, `[-180, 180]`)
#'    - counterclockwise is positive
#'    - clockwise is negative
#' * azimuth angle (in deg, `[-180, 180]`)
#'    - the angle with north (`c(0, 1)`)
#'
#' @examples
#' cal_angle(c(1, 0), c(2, 0)) # == 0
#' cal_angle(c(0, 1), c(0, 2)) # == 0
#'
#' p1 <- c(1, 0)
#' p2 <- c(0, 1)
#'
#' cal_angle(p1, p2) # == 90, counterclockwise
#' cal_angle(p2, p1) # == -90, clockwise
#' @export
cal_angle <- function(v1, v2, clockwise = -1) {
  if (!(is_coord(v1) & is_coord(v2))) stop("v's length or ncol should be 2!")
  # cos = a b / abs(a) * abs(b)
  tmp <- {
    .dot(v1, v2) / (.norm(v1) * .norm(v2))
  } %>% as.numeric()
  sign <- ifelse(sel_x(v2) < sel_x(v1), 1, -1) * (-clockwise)

  # 能否判断是顺时针，还是逆时针？
  rad2deg(acos(tmp) * sign)
}

# 方位角：与正北之间的夹角
# #references
# https://zh.m.wikipedia.org/zh-tw/%E6%96%B9%E4%BD%8D%E8%A7%92

#' @rdname cal_angle
#' @export
cal_azimuth <- function(v, refer = c(0, 1)) {
  ans = cal_angle(refer, v, clockwise = 1)
  ans[ans < 0] %<>% add(360)
  ans
}
