## Micro-DSL (v2.6) & Output Format
Your output must be **Pure Markdown** that implies the structure defined by the Micro-DSL. **Do NOT output EBNF token names** (e.g., `H1_TOKEN`). Generate the actual Markdown (e.g., `# My Title`, `@f my_func(...)`).
**1. Markup Tokens (Implicitly Handled by Markdown):**
`H1` = `# text`
`H2` = `## text`
`NL` = Newline (use sparingly, primarily to end logical entries or separate blocks)
`HR` = `---` (use to separate major logical groups within a section, or between sections)
`Bul` = `- ` (dash + space for general bullets in Header/Legend/Deps)
`IndBul`= `  - ` (two spaces + dash + space for indented `Desc` lines under an `@sigil` entry)
**2. DSL Sigils:**
`@f` = Function `@d` = Data object (e.g., from `data()`)
`@g` = S4/S3 Generic `@m` = S4/S3 Method (class dispatch)
`@c` = R6/R7 constructor (if present)
**3. Other DSL Symbols:**
`|` = Separates multiple function names (e.g., `name1|name2`)
`[...]` = Used for:
* Constructor variants: `ConstructorName[VariantA|VariantB]`
* Method dispatch classes: `methodName[ClassA|ClassB]`
`(...)` = Parameter list in an entry signature.
`param?` = Optional parameter.
`param?=val`= Optional parameter with a default value.
`|` = Separates signature from short description (after params or name if no params).
`->` = Separates short description from return type. Omit if no return value (side-effect).
`!` = Prefix for inline notes (e.g., `!dep`, `!note: text`). There is no space between `!` and the note keyword.
**4. Document Skeleton:**
`Legend?` (H2 "Legend:" + type abbreviation table)
`Header` (H1 PackageName; optional H2 "Core Purpose:", H2 "Key Objects & Concepts:")
`Sections+` (H2 `1.` Title; H2 Title (unnumbered ok); optional H3 "Usage:")
`Entries*` (@sigil lines + optional indented bullets)
`Deps?` (H2 "Key R Dependencies:")
**5. Entry Line Structure:**
`@sigil name(|name)*[variant|ClassA|ClassB]? (alias alt1,alt2)? (param1?=val, param2?, ...)? | Short, pithy description -> ReturnTypeAbbr !note_type: Optional note text`
* **Rules for Entry Line:**
* Omit `()` if no parameters.
* Omit `-> ReturnTypeAbbr` if function has no return value (side-effect only).
* Bundle identical signatures using `name1|name2`.
* Use `ConstructorName[VariantA|VariantB]` for constructor subtypes.
* Use `methodName[DispatchClassA|DispatchClassB]` for S4/S3 methods.
* Notes (`!notetype: text` or `!notetype`) are optional postfixes. Ensure no leading space (e.g., `!ok` not `! ok`).
* Truncate parameter list with `...` if it exceeds 8 parameters (e.g., `(param1, param2, ..., param8, ...)`).
* Example of grouping aliases and optional params: `@f read_csv|read_tsv (file, col_types?="auto", ...) | Parse delimited file -> tib`
**6. Indented Description Bullets (`Desc` lines under an Entry):**
* Format: ` - param_name : type_abbr (constants: val1, "val2" | key_funcs: fnA, fnB)? Brief, essential clarification.`
* Include `(constants: ...)` for params that take a small, fixed set of string literals.
* Include `(key_funcs: ...)` for params that expect specific functions from the package as input.
* **Only include if adding significant clarity** beyond the signature. Omit for common/obvious parameters (e.g., `x`, `...`) or standard defaults (e.g., `drop=TRUE`).
* If an `Entry` already fits in ≤ 110 characters, do not add `Desc` lines unless they would prevent ambiguity.
* Can also be plain text for general notes: ` - General descriptive point.`
**7. Type Abbreviations:**
* Use very short (1-4 letter) type abbreviations (e.g., `NS` for NeuroSpace, `iv` for integer vector, `chr` for character, `log` for logical, `mat` for matrix, `obj` for generic S4/R6 object).
* Reuse type abbreviations whenever identical across entries; do not invent synonyms (e.g., use `int` consistently, not `int`, `intv`, `iv` interchangeably).
**7. Type Abbreviations (deprecated – see next subsection):**
* Use very short (1-4 char) codes, reusing consistently (e.g., `int` not `intv`).
* Built-in codes: int, dbl, num, chr, lgl, lst, vec, df, tib, tbl, mat, arr, fn, env, obj, NS
* Provide a `## Legend:` block only when introducing abbreviations beyond this list.
## Compression Heuristics & Content Selection:
1. **Focus:** Public, user-facing API. Omit internal helpers, unexported symbols, and direct S4 slot accessors (like `slotNames` or methods that just return a slot value if there's already a clear getter). Include only symbols present in NAMESPACE export list.
2. **Grouping:**
* Group trivial getters or functions with identical signatures and purpose using `name1|name2`.
* Group constructors with identical fields but different return subtypes using `ConstructorName[VariantA|VariantB]`.
* Group S3/S4 methods with identical implementations/docs using `methodName[ClassA|ClassB]`.
3. **Methods:** Define generics with `@g`. Emit methods (`@m`) **only** when their behavior, parameters, or return type significantly differ from the generic, or to explicitly list key supported classes.
4. **Omissions:** Skip indented parameter descriptions (`Desc` lines) for obvious defaults (e.g., `drop=TRUE`, `smooth=FALSE`) or very common arguments like `x` or `...` unless they have package-specific meaning. The cheatsheet is not full documentation.
5. **Notes:** Use `!` notes sparingly for critical info (e.g., `!dep`, `!imp`, `!retlist`, `!side`).
6. **Re-exports:** Skip functions and generics re-exported from other packages (e.g., if `dplyr::filter` is re-exported, do not list it).
## Output Contract:
* **Pure Markdown only.** Adhere strictly to the Micro-DSL v2.6.
* No commentary, no intro/outro paragraphs, no code fences unless part of a `CODE_BLOCK_TOKEN` within a `BlockContent` (rarely needed for cheatsheets).
* If any line violates the DSL, regenerate until fully compliant—no prose explanations.
* Generation stops at first line that begins with a second `#` at H1 depth (e.g., if `\\n#` is used as a stop sequence).
* Use a single blank line to separate `Entry` blocks if it aids readability, but avoid excessive blank lines. Use `---` (`HR_TOKEN`) to separate major thematic groups within a section or at the end of sections.
* All parens/brackets/pipes must be balanced.
## Self-Check (Mental Step - Crucial):
Before finalizing, review your output against these critical checks:
1. Does every content line belong to a defined DSL structure (Header, Legend, Section, Entry, Desc, Deps, Block bullet)?
2. Is the DSL syntax for `Entry` lines (sigils, names, params, `|`, `->`, `!`) correctly used? No bare `->` tokens.
3. (reserved)
4. For any abbreviation NOT in the built-in list, is it defined in `## Legend:`?
5. Have you omitted non-essential details and internal functions?

## Few-shot Exemplar
```markdown
# dummyPkg
## Legend:
- int : integer
- chr : character
## 1. Core Functions
@f add (x, y) | Sum two ints -> int
```

## Formal Grammar "v2.6 Micro-EBNF"
Cheatsheet ::= Header Legend? Section+ Deps?
Legend ::= H2_TOKEN TEXT_CONTENT NEWLINE_TOKEN Block
Header ::= H1_TOKEN TEXT_CONTENT NEWLINE_TOKEN+ (H2Section)*
H2Section ::= H2_TOKEN TEXT_CONTENT NEWLINE_TOKEN Block
Section ::= H2_TOKEN (NUMBER_TOKEN PERIOD_TOKEN)? TEXT_CONTENT NEWLINE_TOKEN UsageBlock? Entry+ HR_TOKEN?
UsageBlock ::= H3_TOKEN "Usage:" NEWLINE_TOKEN Block
Entry ::= Sigil_TOKEN EntryIdent ParamList? Bar_TOKEN TEXT_CONTENT ArrowReturn Note? NEWLINE_TOKEN Desc*
Sigil_TOKEN ::= AT_F_TOKEN | AT_D_TOKEN | AT_G_TOKEN | AT_M_TOKEN // Lexer provides @f, @d, @g, @m
EntryIdent ::= IdentGroup MethodOrVariantClass? AliasSpec?
IdentGroup ::= IDENT_TOKEN ("|" IDENT_TOKEN)* // For "foo|bar"
MethodOrVariantClass ::= LBRACKET_TOKEN IDENT_TOKEN (PIPE_TOKEN IDENT_TOKEN)* RBRACKET_TOKEN // For "[ClassA|ClassB]" or "[variantA|variantB]"
AliasSpec ::= LPAREN_TOKEN ALIAS_KEYWORD_TOKEN IDENT_TOKEN (COMMA_TOKEN IDENT_TOKEN)* RPAREN_TOKEN // For "(alias alt1, alt2)"
ParamList ::= LPAREN_TOKEN Param (COMMA_TOKEN Param)* RPAREN_TOKEN
Param ::= IDENT_TOKEN (EQUALS_TOKEN DefaultValue)? OPTIONAL_MARKER_TOKEN?
DefaultValue ::= LITERAL_TOKEN | IDENT_TOKEN
ArrowReturn ::= (ARROW_TOKEN IDENT_TOKEN)?
Note ::= EXCLAMATION_TOKEN NOTETYPE_TOKEN (COLON_TOKEN TEXT_CONTENT)? NEWLINE_TOKEN?
Desc ::= INDENT_TOKEN BULLET_MARKER_TOKEN (ParamDesc | TEXT_CONTENT) NEWLINE_TOKEN
ParamDesc ::= IDENT_TOKEN COLON_TOKEN TYPE_ABBR_TOKEN ParamExtra? TEXT_CONTENT?
ParamExtra ::= LPAREN_TOKEN (ConstantsSpec | KeyFuncsSpec) RPAREN_TOKEN
ConstantsSpec ::= "constants:" (IDENT_TOKEN|LITERAL_TOKEN) (COMMA_TOKEN (IDENT_TOKEN|LITERAL_TOKEN))*
KeyFuncsSpec ::= "key_funcs:" IDENT_TOKEN (COMMA_TOKEN IDENT_TOKEN)*
Deps ::= H2_TOKEN "Key R Dependencies:" NEWLINE_TOKEN Block
Block ::= (Bullet | TEXT_CONTENT | CODE_BLOCK_TOKEN)* (NEWLINE_TOKEN | EOF_TOKEN)
Bullet ::= BULLET_MARKER_TOKEN TEXT_CONTENT NEWLINE_TOKEN
/* --- LEXER-IMPLIED TOKENS (Illustrative) ---
All previous tokens from v2.4, plus:
AT_G_TOKEN, AT_M_TOKEN // @g, @m
// The lexer provides LBRACKET_TOKEN, IDENT_TOKEN, PIPE_TOKEN, RBRACKET_TOKEN.
// The parser, guided by the Sigil_TOKEN (@f for constructor variants, @m for method classes),
// will interpret the content of MethodOrVariantClass appropriately.
*/

---

# neuroim2

## 1. Anatomical Axes & Orientation Constants
### Usage:
- Use these constants to specify anatomical axes or orientations for spatial operations, reorientation, or axis set construction.
- Pass as parameters to functions like `reorient()`, `findAnatomy3D()`, or when constructing axis sets.
@d LEFT_RIGHT | NamedAxis: Left-to-Right anatomical axis -> obj
@d RIGHT_LEFT | NamedAxis: Right-to-Left anatomical axis -> obj
@d ANT_POST | NamedAxis: Anterior-to-Posterior anatomical axis -> obj
@d POST_ANT | NamedAxis: Posterior-to-Anterior anatomical axis -> obj
@d INF_SUP | NamedAxis: Inferior-to-Superior anatomical axis -> obj
@d SUP_INF | NamedAxis: Superior-to-Inferior anatomical axis -> obj
@d TIME | NamedAxis: Time axis for temporal dimension -> obj
@d None | NamedAxis: Null axis (no direction) -> obj
@d NullAxis | AxisSet: Null axis set (no axes) -> obj
@d TimeAxis | AxisSet1D: Time axis set (1D) -> obj
@d OrientationList2D | List of standard 2D anatomical orientations (AxisSet2D) -> lst
@d OrientationList3D | List of standard 3D anatomical orientations (AxisSet3D) -> lst

---

## 2. File Format Constants
### Usage:
- Use these pre-defined objects to specify or detect neuroimaging file formats for reading/writing images and metadata.
- Pass as `descriptor` to meta info constructors or use in file format checks.
@d NIFTI | Standard NIfTI single-file format descriptor -> obj
@d NIFTI_GZ | NIfTI single-file, gzip-compressed format descriptor -> obj
@d NIFTI_PAIR | NIfTI pair (hdr/img) format descriptor -> obj
@d NIFTI_PAIR_GZ | NIfTI pair, gzip-compressed format descriptor -> obj
@d AFNI | AFNI format descriptor -> obj
@d AFNI_GZ | AFNI, gzip-compressed format descriptor -> obj

---

## 3. Core Data Structures & Constructors
### Usage:
- Construct core neuroimaging objects for 2D/3D/4D/5D data using these functions and classes.
- For sparse or masked data, supply a mask (logical array or `LogicalNeuroVol`).
- For file-backed or memory-mapped data, use the corresponding constructor with a file path.
@f NeuroSpace (dim, spacing?=NULL, origin?=NULL, axes?=NULL, trans?=NULL) | Create spatial reference object -> NS
@f NeuroVol (data, space, label?="", indices?=NULL) | Create 3D dense volume -> obj
@f DenseNeuroVol (data, space, label?="", indices?=NULL) | Create dense 3D volume -> obj
@f SparseNeuroVol (data, space, indices?=NULL, label?="") | Create sparse 3D volume -> obj
@f LogicalNeuroVol (data, space, label?="", indices?=NULL) | Create logical (mask) 3D volume -> obj
@f ClusteredNeuroVol (mask, clusters, label_map?=NULL, label?="") | Create clustered 3D volume -> obj
@f NeuroVec (data, space?=NULL, mask?=NULL, label?="") | Create 4D vector (dense or sparse) -> obj
@f DenseNeuroVec (data, space, label?="none") | Create dense 4D vector -> obj
@f SparseNeuroVec (data, space, mask, label?="") | Create sparse 4D vector -> obj
@f BigNeuroVec (data, space, mask, label?="", type?="double", backingfile?=tempfile()) | Create disk-backed sparse 4D vector -> obj
@f FileBackedNeuroVec (file_name, label?=basename(file_name)) | Create file-backed 4D vector -> obj
@f MappedNeuroVec (file_name, label?=basename(file_name)) | Create memory-mapped 4D vector -> obj
@f NeuroVecSeq (...) | Create sequence of NeuroVec objects -> obj
@f NeuroSlice (data, space, indices?=NULL) | Create 2D slice object -> obj
@f IndexLookupVol (space, indices) | Create index lookup volume for sparse mapping -> obj
@f ROICoords (coords) | Create ROI from coordinate matrix -> obj
@f ROIVol (space, coords, data) | Create ROI volume with values -> obj
@f ROIVec (vspace, coords, data?=rep(nrow(coords),1)) | Create ROI vector (4D) -> obj
@f Kernel (kerndim, vdim, FUN?=dnorm, ...) | Create image kernel object -> obj
@f MetaInfo (Dim, spacing, origin?=rep(0,length(spacing)), data_type?="FLOAT", label?="", spatial_axes?=OrientationList3D$AXIAL_LPI, additional_axes?=NullAxis) | Create image meta info -> obj
@f NIFTIMetaInfo (descriptor, nifti_header) | Create NIfTI meta info object -> obj
@f AFNIMetaInfo (descriptor, afni_header) | Create AFNI meta info object -> obj
@f MappedNeuroVecSource (file_name) | Create mapped neurovec source -> obj
@f NeuroVolSource (input, index?=1) | Create source for single 3D volume -> obj
@f NeuroVecSource (file_name, indices?=NULL, mask?=NULL) | Create source for 4D vector -> obj
@f SparseNeuroVecSource (meta_info, indices?=NULL, mask) | Create source for sparse 4D vector -> obj

---

## 4. Pre-defined Orientation & Axis Sets
### Usage:
- Use these lists to select standard anatomical orientations for constructing or reorienting spaces.
- Pass elements from `OrientationList2D` or `OrientationList3D` as axes or orientation parameters.
@d OrientationList2D | List of standard 2D anatomical orientations (AxisSet2D) -> lst
@d OrientationList3D | List of standard 3D anatomical orientations (AxisSet3D) -> lst

---

## 5. File Format & Metadata Utilities
### Usage:
- Use these functions to check file format compatibility, extract or strip file extensions, and read meta info.
@f file_matches (x, file_name) | Check if file matches format -> lgl
@f header_file_matches (x, file_name) | Check if file is header for format -> lgl
@f data_file_matches (x, file_name) | Check if file is data for format -> lgl
@f header_file (x, file_name) | Get header file name for format -> chr
@f data_file (x, file_name) | Get data file name for format -> chr
@f strip_extension (x, file_name) | Remove extension from file name -> chr
@f read_meta_info (x, file_name) | Read meta info for file/format -> obj
@f read_header (file_name) | Read image header, auto-detect format -> obj

---

## 6. ROI & Patch Constructors
### Usage:
- Use these helpers to create ROIs of various shapes for analysis or searchlight methods.
@f square_roi (bvol, centroid, surround, fill?=NULL, nonzero?=FALSE, fixdim?=3) | Create square ROI at centroid -> obj
@f cuboid_roi (bvol, centroid, surround, fill?=NULL, nonzero?=FALSE) | Create cuboid ROI at centroid -> obj
@f spherical_roi (bvol, centroid, radius, fill?=NULL, nonzero?=FALSE, use_cpp?=TRUE) | Create spherical ROI at centroid -> obj
@f spherical_roi_set (bvol, centroids, radius, fill?=NULL, nonzero?=FALSE) | Create multiple spherical ROIs -> lst

---

## 7. Searchlight & Clustering Utilities
### Usage:
- Use these to generate searchlight iterators or cluster-based ROI sets for multivariate analyses.
@f random_searchlight (mask, radius) | Random searchlight iterator (spherical ROIs) -> lst
@f bootstrap_searchlight (mask, radius?=8, iter?=100) | Bootstrap searchlight iterator -> lst
@f searchlight_coords (mask, radius, nonzero?=FALSE, cores?=0) | Exhaustive searchlight (coords) -> lst
@f searchlight (mask, radius, eager?=FALSE, nonzero?=FALSE, cores?=0) | Exhaustive searchlight (ROIs) -> lst
@f clustered_searchlight (mask, cvol?=NULL, csize?=NULL) | Clustered searchlight iterator -> lst

---

## 8. Spatial Filtering & Enhancement
### Usage:
- Apply spatial smoothing or edge-preserving filters to volumetric data.
@f gaussian_blur (vol, mask, sigma?=2, window?=1) | Gaussian blur on volume -> obj
@f guided_filter (vol, radius?=4, epsilon?=0.49) | Edge-preserving guided filter -> obj
@f bilateral_filter (vol, mask, spatial_sigma?=2, intensity_sigma?=1, window?=1) | Bilateral filter on volume -> obj
@f laplace_enhance (vol, mask, k?=2, patch_size?=3, search_radius?=2, h?=0.7, mapping_params?=NULL, use_normalization_free?=TRUE) | Laplacian enhancement filter -> obj

---

## 9. Binary I/O Utilities
### Usage:
- Use these to create binary readers/writers for custom data access.
@f BinaryReader (input, byte_offset, data_type, bytes_per_element, endian?=.Platform$endian, signed?=TRUE) | Create binary data reader -> obj
@f BinaryWriter (output, byte_offset, data_type, bytes_per_element, endian?=.Platform$endian) | Create binary data writer -> obj
@f ColumnReader (nrow, ncol, reader) | Create column-wise data reader -> obj

---

## 10. Miscellaneous
### Usage:
- Use these for quaternion/affine conversions or kernel embedding.
@f matrixToQuatern (mat) | Convert affine matrix to quaternion -> lst
@f quaternToMatrix (quat, origin, stepSize, qfac) | Convert quaternion to affine matrix -> mat
@f embed_kernel (x, sp, center_voxel, ...) | Embed kernel in space at voxel -> obj

---

## 11. Data Object Families (Constants for Orientation)
### Usage:
- Use these constants as choices for orientation parameters or when constructing axis sets.
@d LEFT_RIGHT | NamedAxis: Left-to-Right anatomical axis -> obj
@d RIGHT_LEFT | NamedAxis: Right-to-Left anatomical axis -> obj
@d ANT_POST | NamedAxis: Anterior-to-Posterior anatomical axis -> obj
@d POST_ANT | NamedAxis: Posterior-to-Anterior anatomical axis -> obj
@d INF_SUP | NamedAxis: Inferior-to-Superior anatomical axis -> obj
@d SUP_INF | NamedAxis: Superior-to-Inferior anatomical axis -> obj

---

## 12. Null/None Axis Constants
### Usage:
- Use these as placeholders for missing or undefined axes.
@d None | NamedAxis: Null axis (no direction) -> obj
@d NullAxis | AxisSet: Null axis set (no axes) -> obj

---

## 13. Time Axis Constants
### Usage:
- Use these for temporal axes in axis sets or NeuroSpace.
@d TIME | NamedAxis: Time axis for temporal dimension -> obj
@d TimeAxis | AxisSet1D: Time axis set (1D) -> obj

---

## 14. Pre-defined Orientation Lists
### Usage:
- Use these lists to select standard anatomical orientations for constructing or reorienting spaces.
@d OrientationList2D | List of standard 2D anatomical orientations (AxisSet2D) -> lst
@d OrientationList3D | List of standard 3D anatomical orientations (AxisSet3D) -> lst

---

## 15. File Format Descriptor Constants
### Usage:
- Use these constants to specify or detect file formats for reading/writing images and metadata.
@d NIFTI | Standard NIfTI single-file format descriptor -> obj
@d NIFTI_GZ | NIfTI single-file, gzip-compressed format descriptor -> obj
@d NIFTI_PAIR | NIfTI pair (hdr/img) format descriptor -> obj
@d NIFTI_PAIR_GZ | NIfTI pair, gzip-compressed format descriptor -> obj
@d AFNI | AFNI format descriptor -> obj
@d AFNI_GZ | AFNI, gzip-compressed format descriptor -> obj

---

#