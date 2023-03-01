

# This is a function to check if an object is an OTS

# Input parameters
# X: an object

#--------------------------------------------------------------------------------------

check_ots <- function(X) {


  # The element must be a numeric vector


  if (!is.numeric(X)) {

  stop('The object must be a numeric vector')

  }


  # The element can not contain NA entries

  check_nas <- sum(is.na(X))

  if (sum(check_nas) != 0) {

    stop('There are some NAs in the series')

  }

}
