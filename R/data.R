#' Fonville et al. haemagglutination inhibition titer data
#'
#' Longitudinal HI titer measurements from Fonville et al. (2014).
#' Participants were sampled across multiple years; use
#' \code{group_col = "Sample Year"} when computing per-sample metrics.
#' Data are in long format with one row per participant-strain-sample combination.
#'
#' @format A data frame with 9 columns:
#' \describe{
#'   \item{Subject Number}{Participant identifier.}
#'   \item{Sample Year}{Year the serum sample was collected (numeric).}
#'   \item{Year of Birth}{Participant's year of birth (numeric).}
#'   \item{Sample}{Serum sample identifier.}
#'   \item{PCR Results}{PCR test result recorded for the sample.}
#'   \item{Row in Fonville Fig S15}{Row index in Figure S15 of the source
#'     publication.}
#'   \item{strain}{Influenza strain name (character).}
#'   \item{titer}{Haemagglutination inhibition titer (numeric). Below-detection
#'     and censored values are recoded during data preparation.}
#'   \item{isolation_year}{Year the test strain was isolated (numeric).}
#' }
#'
#' @source Fonville et al. (2014) Antibody landscapes after influenza virus
#'   infection or vaccination. \emph{Science}, 346(6212), 996-1000.
#'   \doi{10.1126/science.1256427}
#'
#' @examples
#' data(Fonville)
#' head(Fonville)
"Fonville"


#' Lessler et al. haemagglutination inhibition titer data
#'
#' HI titer measurements from Lessler et al. (2012). Each participant
#' (\code{id}) has titers against a panel of influenza strains. Several
#' strains share the same isolation year, so use \code{mode = "mean"}
#' (or "max"/"min") to aggregate duplicate isolation years when computing
#' distance-weighted metrics. Data are in long format with one row per
#' participant-strain combination.
#'
#' @format A data frame with 9 columns:
#' \describe{
#'   \item{age}{Participant age in years (numeric).}
#'   \item{is.vac}{Vaccination status indicator.}
#'   \item{shift.age}{Age-related covariate from the source dataset (numeric).}
#'   \item{titers}{Titer summary column from the source dataset.}
#'   \item{neut.against}{Strain against which neutralisation was measured.}
#'   \item{id}{Participant identifier.}
#'   \item{strain}{Influenza strain name (character).}
#'   \item{titer}{Haemagglutination inhibition titer (numeric).}
#'   \item{isolation_year}{Year the test strain was isolated (numeric).}
#' }
#'
#' @source Lessler et al. (2012) Evidence for antigenic seniority in influenza
#'   A (H3N2) antibody responses in southern China. \emph{PLOS Pathogens},
#'   8(7), e1002802. \doi{10.1371/journal.ppat.1002802}
#'
#' @examples
#' data(Lessler)
#' head(Lessler)
"Lessler"
