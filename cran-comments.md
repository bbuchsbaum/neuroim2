## R CMD check results

* 0 errors ✔ | 0 warnings ✔ | 0 notes ✖

* This is a new release.

* We have added a reference in description

* We have removed "dontrun" from examples and instead used "donttest" in some rare cases

* updated spelling wordlist to avoid mispelling NOTE

* We have reset user options that were set with "par".

* We have made the title shorter than 65 characters.

* devtools::check_win_devel() reveals mispelled words in description, but these are false alarms.

* removed acronyms in title description ("fMRI" is used only after earlier definition)

* After hunting for "print/cat" statements in non-display code, we found 1 and deleted it. 

* Many print statements remain in S4 "show" methods and in @examples. We have left these in place.

* We have added @return tags for all remaining functions that we could find.

