## Micro-DSL v3.4 Source Format Grammar

This is the HUMAN-READABLE source format for v3.4. It extends v3.3 with semantic annotations
and constrained types while maintaining clarity and completeness. A separate compilation step 
produces the compressed format.

**Budget**: Target ≤ 250 lines or ≤ 2,500 tokens for optimal LLM processing.

**1. Document Structure:**
```
@pkg package_name | description
[Type Aliases section]
[Constraint Definitions section]
[Legend section if needed]
# Package Name
[Sections with entries]
[Dependencies section]
[Meta-Footer]
```

**2. Sigils (Same as v3.3):**
```
@pkg - Package declaration
@f   - Function
@d   - Data object
@x   - Re-export from another package

S3 System:
@s3g - S3 generic (UseMethod)
@s3m - S3 method (generic.class)  
@s3c - S3 class definition

S4 System:
@s4g - S4 generic (setGeneric)
@s4m - S4 method (setMethod)
@s4c - S4 class (setClass)

S7 System:
@s7g - S7 generic (new_generic)
@s7m - S7 method
@s7c - S7 class (new_class)

R6 System:
@r6c - R6 class (R6Class)
```

**3. Type System (v3.4 Enhanced with Constraints):**
```
Scalars (default): int, dbl, chr, lgl, raw, cpl
Vectors: vec<type> or type[] 
Matrices: mat<type> or mat<type,rows,cols>
Arrays: arr<type,dims>
Lists: lst<type> or lst{field:type, ...} (structured)
Data frames: df, tbl, data.table
Factors: fct, ord
Dates: Date, POSIXct, POSIXlt

Union types: type1|type2|type3
Nullable: type? (shorthand for type|NULL)
Any type: any
Ellipsis: ... or ...:type (e.g., ...:expr for NSE)

Class types: s3:classname, s4:classname, r6:classname, s7:classname

CONSTRAINED TYPES (v3.4):
Enums: chr["opt1"|"opt2"|"opt3"]
Ranges: int[min..max], dbl[min..max]
Patterns: chr[/regex/]
Exclusions: int[1..100]&!=[13,17]
References: @ref:constraint_name
```

**4. Entry Format (v3.4 Enhanced):**
```
@sigil name (param1:type1[constraint]?=default, param2:type2, ...) | Description -> return_type | return_schema
  +tag:value +tag:value
  !cov [Class1, Class2] (for generics)
  - param1 : Additional description
    @annotation:value @annotation:value
  - param2 : (constants: "a", "b", CONST) Valid values
    @requires:condition @affects:target
  - param3 : (key_funcs: func1, func2) Related functions
    @lifecycle:init @units:measurement
  ```verbatim
  # Optional verbatim R code block
  ```
```

**5. Type Aliases Section (v3.4 Enhanced):**
```markdown
## Type Aliases:
DF = df|tbl|data.table              # Standard data frame types
V<T> = vec<T>                       # Vector shorthand
Fml = s3:formula                    # Formula objects
Gg = s3:ggplot                      # ggplot2 objects
Config = lst{method:chr, opts:lst}  # Structured config
ValidPort = int[1024..65535]       # Constrained port range
```

**Standard aliases** (use these by default):
- `DF` for data frame arguments
- `Fml` for formula arguments (not `fml`)
- `V<T>` for vectors when brevity helps

**6. Constraint Definitions (v3.4 New):**
```markdown
## Constraint Definitions:
@constraint positive_weights | Positive numeric weights
  type: vec<dbl>
  validates: all(. > 0)
  length: @env:nrow(data)
  
@constraint valid_identifier | Valid R identifier
  type: chr
  pattern: /^[a-zA-Z_][a-zA-Z0-9_.]*$/
  not_reserved: TRUE
```

**7. Class Documentation (v3.4 Enhanced):**
```
@s4c ClassName | One-line description
  - slots: name1 (type1[constraint]) @annotation:value
           name2 (type2) @lifecycle:init @immutable
  - extends: ParentClass
  - validity: Description of validity rules
  
@r6c ClassName | One-line description  
  - fields: field1 (type1[constraint]) @purpose:role
            field2 (type2) @lazy @cached
  - methods: method1 (args) -> ret_type
             method2 (args) -> ret_type  
  - inherits: ParentClass
```

**8. Metadata Tags (v3.3 + v3.4 additions):**
```
+family:group_name           # Function family
+pipe:in|out                # Pipe compatibility (in, out, or both)
+nse:param1,param2          # Parameters using NSE
+side:effect[details]       # Side effects with sub-facets
  - fs[read|write|delete]   # File system operations
  - plot[device|file]       # Graphics output
  - console[print|message|warning]  # Console output
  - network[http|socket|download]   # Network operations
  - options[get|set|env]    # Global options/environment
  - db[read|write|query]    # Database operations
+perf:O(complexity)         # Performance complexity
+mem:usage                  # Memory usage pattern
+compute:intensity          # Computational intensity
+deprecated:replacement     # Deprecation with suggested alternative
+wraps:function            # This function wraps another
+calls:func1,func2         # Functions called internally
+see:related1,related2     # Related functions to consider
+parallel:capable          # Can use parallel processing (v3.4)
+deterministic:false       # Non-deterministic results (v3.4)
+pure:false               # Has side effects (v3.4)
```

**9. Semantic Annotations (v3.4 New):**
```
BEHAVIORAL:
@controls:aspect          # Parameter controls specific behavior
@affects:target          # Changes affect another component  
@modifies:target         # Directly modifies target

DEPENDENCY:
@requires:condition      # Prerequisite condition
@conflicts:parameter     # Mutually exclusive with
@extends:base           # Extends functionality

VALIDATION:
@validates:constraint    # Validation rule
@range:[min,max]        # Numeric range
@length:constraint      # Length requirement
@pattern:regex          # Pattern matching

SEMANTIC ROLE:
@purpose:role           # Semantic purpose
@units:measurement      # Physical/logical units
@example:value          # Example values
@default-reason:why     # Why this default

LIFECYCLE:
@lifecycle:stage        # When relevant (init|config|runtime|cleanup)
@immutable             # Cannot be modified
@cached                # Result is cached
@lazy                  # Evaluated on demand

CONDITIONAL:
@when:condition        # Conditional applicability
@implies:consequence   # Logical implication
@if:cond @then:result  # If-then constraints
```

**10. Structured Return Types (v3.4 New):**
```
-> lst{
  field1: type1 @annotation,
  field2: type2[constraint] @annotation,
  nested: lst{
    subfield: type
  }
}
```

**11. Example Entry (v3.4):**
```markdown
@f analyze_model (
  model:s3:lm,
  type:chr["summary"|"anova"|"diagnostics"]?="summary",
  conf.level:dbl[0.5..0.99]?=0.95
) | Analyze fitted model -> lst{
  statistics: df @purpose:results,
  plots: lst<s3:ggplot>? @when:type="diagnostics",
  interpretation: chr @purpose:summary
}
  +family:analysis +compute:light
  - model : @requires:fitted @validates:has-residuals
  - type : @controls:output-format @affects:return-structure
  - conf.level : @purpose:confidence @affects:statistics.ci
```

**12. Conditional Constraints (v3.4):**
```markdown
@f process_data (
  data:df,
  method:chr["scale"|"center"|"none"]?="none",
  scale.center:lgl?=TRUE,
  scale.scale:lgl?=TRUE
) | Process data with scaling options -> df
  - method : @controls:processing
  - scale.center : @when:method="scale" @requires:TRUE
                   @when:method="center" @implies:scale.scale=FALSE
  - scale.scale : @when:method="scale" @default:TRUE
                  @conflicts:method="center"
```

**13. Best Practices (v3.4):**
- Use specific sigils (@s3g not @g)
- Always specify vector vs scalar types
- Use standard type aliases (DF, Fml, V<T>)
- Add constraints from match.arg/stopifnot/checks
- Keep !cov lists short (3-6 classes max)
- Document semantic relationships concisely
- Use structured types for complex returns
- Define reusable constraints with @constraint
- Include conditional logic with @when/@implies
- Group related functions with +family tags
- Mark side effects with detailed sub-facets
- Stay within budget (≤250 lines)

**14. Meta-Footer:**
```markdown
---
## Meta-Footer
- Micro-DSL Version: v3.4-source
- Package: {pkg} (Version: X.Y.Z)
- Generated: [ISO-8601 timestamp]
- Features: types[constrained] sigils[specific] metadata[rich] semantics[annotated]
- Coverage: {n_documented_exports} / {n_total_exports} exports
- Provenance: exports[NAMESPACE], enums[match.arg/switch], constraints[assertions/checks]
```

**15. Export Detection Priority:**
1. NAMESPACE file: `export()`, `S3method()`, `exportClasses()`, `exportMethods()`, `exportPattern()`
2. Roxygen tags: `@export` in documentation
3. If neither present: skip the symbol (do not guess or include)

**16. Inference Heuristics (apply silently):**
- Type from defaults: TRUE/FALSE → lgl, "text" → chr, 1L → int, 1.0 → dbl
- Common patterns: data/df/tbl → DF, formula → Fml, weights → vec<dbl>
- Enums: match.arg(x, c("a","b")) → chr["a"|"b"]
- Ranges: stopifnot(x >= 0 && x <= 1) → dbl[0..1]
- Side effects: file.* → fs, plot/ggplot → plot, message/cat → console
- Determinism: runif/sample/rnorm → +deterministic:false

---

```markdown
@pkg neuroim2 | Data Structures for Brain Imaging Data

## Type Aliases:
DF = df|tbl|data.table
V<T> = vec<T>
Fml = s3:formula
Gg = s3:ggplot

## Constraint Definitions:
@constraint positive_weights | Positive numeric weights
  type: vec<dbl>
  validates: all(. > 0)

# neuroim2

## 1. Core

@f ANT_POST () | Pre-defined anatomical axis -> s3:NamedAxis
  +family:anatomical_axes

@f BigNeuroVec (data:arr<dbl,4>, space:s4:NeuroSpace, mask:s4:LogicalNeuroVol, label:chr?="", type:chr["double"|"float"|"integer"]?="double", backingfile:chr?=tempfile()) | Create memory-mapped neuroimaging vector -> s4:BigNeuroVec
  +family:neurovec +pipe:in|out +compute:moderate
  - data : @purpose:input-data
  - space : @purpose:spatial-properties
  - mask : @purpose:voxel-selection
  - label : @purpose:label
  - type : @purpose:storage-type
  - backingfile : @purpose:memory-mapping

@f BinaryReader (input:chr|con, byte_offset:int, data_type:chr, bytes_per_element:int, endian:chr?=.Platform$endian, signed:lgl?=TRUE) | Create binary reader -> s4:BinaryReader
  +family:binary_io +pipe:in
  - input : @purpose:input-source
  - byte_offset : @purpose:offset
  - data_type : @purpose:data-type
  - bytes_per_element : @purpose:element-size
  - endian : @purpose:endian
  - signed : @purpose:signedness

@f BinaryWriter (output:chr|con, byte_offset:int, data_type:chr, bytes_per_element:int, endian:chr?=.Platform$endian) | Create binary writer -> s4:BinaryWriter
  +family:binary_io +pipe:out
  - output : @purpose:output-destination
  - byte_offset : @purpose:offset
  - data_type : @purpose:data-type
  - bytes_per_element : @purpose:element-size
  - endian : @purpose:endian

@f ClusteredNeuroVec (x:s4:NeuroVec|mat<dbl>, cvol:s4:ClusteredNeuroVol, FUN:f?=mean, weights:vec<dbl>?, label:chr?="") | Create clustered neuroimaging vector -> s4:ClusteredNeuroVec
  +family:neurovec +pipe:in|out +compute:moderate
  - x : @purpose:input-data
  - cvol : @purpose:cluster-assignments
  - FUN : @purpose:aggregation-function
  - weights : @purpose:voxel-weights
  - label : @purpose:label

@f ColumnReader (nrow:int, ncol:int, reader:f) | Create column reader -> s4:ColumnReader
  +family:binary_io +pipe:in
  - nrow : @purpose:number-of-rows
  - ncol : @purpose:number-of-columns
  - reader : @purpose:reader-function

@f INF_SUP () | Pre-defined anatomical axis -> s3:NamedAxis
  +family:anatomical_axes

@f LEFT_RIGHT () | Pre-defined anatomical axis -> s3:NamedAxis
  +family:anatomical_axes

@f LogicalNeuroVol (data:arr<lgl,3>, space:s4:NeuroSpace) | Create logical neuroimaging volume -> s4:LogicalNeuroVol
  +family:neurovol +pipe:in|out +compute:moderate
  - data : @purpose:input-data
  - space : @purpose:spatial-properties

@f POST_ANT () | Pre-defined anatomical axis -> s3:NamedAxis
  +family:anatomical_axes

@f RIGHT_LEFT () | Pre-defined anatomical axis -> s3:NamedAxis
  +family:anatomical_axes

@f SUP_INF () | Pre-defined anatomical axis -> s3:NamedAxis
  +family:anatomical_axes

@f TIME () | Pre-defined time axis -> s3:NamedAxis
  +family:anatomical_axes

@f TimeAxis () | Pre-defined time axis set -> s3:AxisSet1D
  +family:anatomical_axes

@s4g resample (source:s4:NeuroVol, target:s4:NeuroSpace, ...:expr) | Resample image to match target space -> s4:NeuroVol
  +family:resample +pipe:in|out +compute:moderate +nse:...
  !cov [NeuroVol, NeuroVec]

@s4g scale (x:s4:NeuroVol, ...:expr) | Scale image -> s4:NeuroVol
  +family:scale +pipe:in|out +compute:moderate +nse:...
  !cov [NeuroVol, NeuroVec]

@s4g slice (x:s4:NeuroVol, zlevel:int, along:chr, orientation:chr, ...:expr) | Extract 2D slice -> s4:NeuroSlice
  +family:slice +pipe:in|out +compute:moderate +nse:...
  !cov [NeuroVol, NeuroVec]

@s4g space (x:s4:NeuroVol, ...:expr) | Extract spatial properties -> s4:NeuroSpace
  +family:space +pipe:in +compute:light +nse:...
  !cov [NeuroVol, NeuroVec]

@s4g values (x:s4:NeuroVol, ...:expr) | Extract values -> vec<dbl>
  +family:values +pipe:in +compute:light +nse:...
  !cov [NeuroVol, NeuroVec]

## Dependencies
- Imports: Rcpp, RcppArmadillo, RcppParallel, Matrix, RNifti, dbscan, stringr, methods, bigstatsr, RNiftyReg, future.apply, deflist, crayon, ggplot2, magrittr, grid
- Suggests: testthat, covr, knitr, roxygen2, rmarkdown, Gmedian, R.utils, spelling

---
## Meta-Footer
- Micro-DSL Version: v3.4-source
- Package: neuroim2 (Version: 0.8.2)
- Generated: 2023-10-05T12:00:00Z
- Features: types[constrained] sigils[specific] metadata[rich] semantics[annotated]
- Coverage: 20 / 100 exports
- Provenance: exports[NAMESPACE], enums[match.arg/switch], constraints[assertions/checks]
```