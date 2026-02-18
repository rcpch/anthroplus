#' Compute measurements against z scores for age 5 to 19
#'
#' @param sex A numeric or text variable containing gender information.
#'            If it is numeric, its values must be: 1 for males and 2 for
#'            females. If it is character, it must be "m" or "M" for males
#'            and "f" or "F" for females. No z-scores will be calculated
#'            if sex is missing.
#' @param age_in_months A numeric variable containing age information;
#'            Age-related z-scores will NOT be calculated if age is missing.
#' @param measurement_method Character: one of `"length"`, `"weight"`, `"bmi"`, `"headc"`.
#' @param requested_z Numeric vector of requested z-scores to invert.
#' @param is_age_in_month Logical; whether `age` is in months (default FALSE).
#' @param measurement_precision Integer digits to round the returned measurement to (default 2L).
#' @param correct_extreme Logical; if TRUE apply the same linear extreme-value correction
#'   used in `compute_zscore_adjusted()` for |z| > 3. Defaults to FALSE.
#' 
#' @return A data.frame with a single column named according to the requested
#' measurement (e.g. `lenhei`, `weight`, `bmi`, `headc`) containing the computed
#' measurement values (or NA where not computable).
#' @examples
#' anthro_measurements(sex = 1, age = 365, requested_z = 0, measurement_method = "length")


#' Invert z-scores to measurements for age 5-19 (anthroplus)
#'
#' @export
anthroplus_measurements <- function(sex,
																		age_in_months = NA_real_,
																		measurement_method = c("length", "weight", "bmi"),
																		requested_z = NA_real_,
																		measurement_precision = 2L,
																		correct_extreme = FALSE) {
	stopifnot(all(tolower(sex) %in% c("1", "2", "f", "m", NA_character_)))
	stopifnot(is.numeric(age_in_months))
	stopifnot(is.numeric(requested_z))
	measurement_method <- match.arg(measurement_method)
	stopifnot(is.logical(correct_extreme) && length(correct_extreme) == 1L)

	n <- max(length(sex), length(age_in_months), length(requested_z))
	sex <- rep_len(sex, n)
	age_in_months <- rep_len(age_in_months, n)
	requested_z <- rep_len(requested_z, n)

	# prepare ages for interpolation like zscore_indicator
	low_age <- trunc(age_in_months)
	upp_age <- trunc(age_in_months + 1)
	diff_age <- age_in_months - low_age
	data <- data.frame(sex = sex, low_age = low_age, upp_age = upp_age, ordering = seq_len(n))

	if (measurement_method %in% c("length", "len", "height")) {
		gs <- hfa_growth_standards
		out_name <- "height_in_cm"
		age_upper_bound <- 229
		age_lower_bound <- 60
	} else if (measurement_method %in% c("weight", "wei")) {
		gs <- wfa_growth_standards
		out_name <- "weight_in_kg"
		age_upper_bound <- 121
		age_lower_bound <- 60
	} else if (measurement_method == "bmi") {
		gs <- bfa_growth_standards
		out_name <- "bmi"
		age_upper_bound <- 229
		age_lower_bound <- 60
	} else {
		stop("Unsupported measurement_method")
	}

	match_low_age <- merge(data, gs, by.x = c("sex", "low_age"), by.y = c("sex", "age"), all.x = TRUE, sort = FALSE)
	match_upp_age <- merge(data, gs, by.x = c("sex", "upp_age"), by.y = c("sex", "age"), all.x = TRUE, sort = FALSE)
	match_low_age <- match_low_age[order(match_low_age$ordering), ]
	match_upp_age <- match_upp_age[order(match_upp_age$ordering), ]

	m <- match_low_age[["m"]]
	l <- match_low_age[["l"]]
	s <- match_low_age[["s"]]

	is_diff_age_pos <- !is.na(diff_age) & diff_age > 0
	if (any(is_diff_age_pos)) {
		adjust_param <- function(x) {
			x_name <- as.character(substitute(x))
			x[is_diff_age_pos] + diff_age[is_diff_age_pos] * (match_upp_age[[x_name]][is_diff_age_pos] - x[is_diff_age_pos])
		}
		m[is_diff_age_pos] <- adjust_param(m)
		l[is_diff_age_pos] <- adjust_param(l)
		s[is_diff_age_pos] <- adjust_param(s)
	}

	# inverse LMS
	calc_sd <- function(mv, lv, sv, sd) {
		ifelse(lv == 0, mv * exp(sv * sd), mv * ((1 + lv * sv * sd)^(1 / lv)))
	}

	y <- rep(NA_real_, n)
	valid_idx <- !is.na(m) & !is.na(l) & !is.na(s) & !is.na(age_in_months)

	# non-extreme
	non_extreme <- valid_idx & (requested_z <= 3 & requested_z >= -3)
	if (any(non_extreme)) {
		lv <- l[non_extreme]
		mv <- m[non_extreme]
		sv <- s[non_extreme]
		zv <- requested_z[non_extreme]
		y[non_extreme] <- ifelse(lv == 0, mv * exp(sv * zv), mv * ((1 + lv * sv * zv)^(1 / lv)))
	}

	if (isTRUE(correct_extreme)) {
		pos_ext <- valid_idx & (requested_z > 3)
		if (any(pos_ext)) {
			mv <- m[pos_ext]
			lv <- l[pos_ext]
			sv <- s[pos_ext]
			SD3pos <- calc_sd(mv, lv, sv, 3)
			SD2pos <- calc_sd(mv, lv, sv, 2)
			SD23pos <- SD3pos - SD2pos
			y[pos_ext] <- SD3pos + (requested_z[pos_ext] - 3) * SD23pos
		}
		neg_ext <- valid_idx & (requested_z < -3)
		if (any(neg_ext)) {
			mv <- m[neg_ext]
			lv <- l[neg_ext]
			sv <- s[neg_ext]
			SD3neg <- calc_sd(mv, lv, sv, -3)
			SD2neg <- calc_sd(mv, lv, sv, -2)
			SD23neg <- SD2neg - SD3neg
			y[neg_ext] <- SD3neg + (requested_z[neg_ext] + 3) * SD23neg
		}
	}

	# enforce age bounds
	valid_age <- !is.na(age_in_months) & age_in_months >= age_lower_bound & age_in_months < age_upper_bound
	y[!valid_age] <- NA_real_

	y <- round(y, digits = as.integer(measurement_precision))
	out <- as.data.frame(matrix(nrow = n, ncol = 0))
	out[[out_name]] <- y
	out
}
