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
#' (\code{id}) was tested once against a panel of 9 H3N2 strains spanning
#' 1968-2008, so each (id, isolation_year) pair appears exactly once and
#' no \code{group_col} or \code{mode} aggregation is needed. Titers are
#' already log-transformed. Data are in long format with one row per
#' participant-strain combination.
#'
#' @format A data frame with 7 columns:
#' \describe{
#'   \item{id}{Participant identifier (numeric).}
#'   \item{age}{Participant age in years (numeric).}
#'   \item{is.vac}{Vaccination status indicator ("1" or NA).}
#'   \item{shift.age}{Age relative to the antigenic-shift year used in the
#'     source publication (numeric).}
#'   \item{strain}{Influenza strain name (character), e.g. "A/HK/1968(H3N2)".}
#'   \item{isolation_year}{Year the test strain was isolated (numeric).}
#'   \item{titer}{Log-transformed haemagglutination inhibition titer (numeric).}
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
