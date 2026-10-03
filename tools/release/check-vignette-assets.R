# Verify the source tarball and installed docs independently. HTML keeps plots,
# scripts and layout CSS embedded; only the unchanged font files are shared.
check_vignette_assets <- function(doc_dir, source_dir = "vignettes") {
  expected <- sub("[.]Rmd$", ".html", list.files(source_dir, "[.]Rmd$"))
  html <- file.path(doc_dir, expected)
  stopifnot(length(expected) > 0L, all(file.exists(html)))
  font_css <- file.path(doc_dir, "albers-fonts.css")
  stopifnot(file.exists(font_css))
  css <- paste(readLines(font_css, warn = FALSE), collapse = "\n")
  font_urls <- regmatches(css, gregexpr('url\\("[^"]+"\\)', css))[[1]]
  font_urls <- sub('^url\\("', "", sub('"\\)$', "", font_urls))
  expected_fonts <- list.files(file.path(source_dir, "fonts"), "[.]woff2$")
  stopifnot(length(font_urls) == 7L,
            setequal(font_urls, file.path("fonts", expected_fonts)))
  resources <- c("albers-fonts.css", font_urls, "fonts/LICENSE")
  installed <- file.path(doc_dir, resources)
  source <- file.path(source_dir, resources)
  stopifnot(all(file.exists(installed)),
            identical(unname(tools::md5sum(installed)),
                      unname(tools::md5sum(source))))
  for (file in html) {
    page <- xml2::read_html(file)
    assets <- xml2::xml_text(xml2::xml_find_all(page, paste(
      "//link/@href | //script/@src | //img/@src | //source/@src |",
      "//object/@data | //iframe/@src")))
    # No network or per-vignette directories are needed to display these docs.
    stopifnot(sum(assets == "albers-fonts.css") == 1L,
              all(assets == "albers-fonts.css" | startsWith(assets, "data:")))
    styles <- paste(xml2::xml_text(xml2::xml_find_all(page, "//style")),
                    collapse = "\n")
    stopifnot(!grepl("@font-face", styles, fixed = TRUE))
    urls <- regmatches(styles, gregexpr("url\\([^)]*\\)", styles))[[1]]
    if (length(urls)) stopifnot(all(grepl("^url\\(['\"]?data:", urls)))
  }
  files <- list.files(doc_dir, recursive = TRUE, full.names = TRUE)
  total <- sum(file.info(files)$size)
  # CRAN's general documentation limit is 5 MB; use decimal MB conservatively.
  stopifnot(total < 5000000)
  data.frame(directory = doc_dir, html_files = length(html), fonts = length(font_urls),
             bytes = total, offline_resources = TRUE, source_assets_identical = TRUE)
}
