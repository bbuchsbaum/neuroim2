# Golden Tests: AI Reference

## Purpose
Golden tests ensure semantic equivalence across language implementations (R, Python, Rust) of the same software by validating numeric outputs rather than implementation details.

## Core Concepts

### What Golden Tests Are
- Language-agnostic test specifications in XML format
- Focus on WHAT code should do (semantic behavior), not HOW
- Validate using numeric outputs with tolerances
- Each language implements same semantics differently

### Key Principles
1. **Numeric Focus**: All validations based on matrix dimensions, values, statistical properties
2. **Semantic Descriptions**: Purpose + mathematical algorithm in each test
3. **Progressive Enhancement**: Spec → R → Python → Rust implementations
4. **Language Agnosticism**: Describe behavior, not implementation

## XML Test Structure

```xml
<?xml version="1.0" encoding="UTF-8"?>
<golden_test xmlns="http://golden-tests.org/schema">
  <metadata>
    <id>unique_test_id</id>
    <version>1.0</version>
    <description>Brief description</description>
    <tags><tag>category</tag></tags>
  </metadata>
  
  <semantic_description>
    <purpose>What functionality is tested</purpose>
    <algorithm>Step-by-step mathematical description</algorithm>
  </semantic_description>
  
  <inputs>
    <!-- Structured test data -->
  </inputs>
  
  <expected_outputs>
    <numeric_checks>
      <check>
        <type>exact_value|range|relative|statistical</type>
        <location>where to check</location>
        <expected>value</expected>
        <tolerance>acceptable deviation</tolerance>
      </check>
    </numeric_checks>
  </expected_outputs>
  
  <implementations>
    <R><![CDATA[# R code]]></R>
    <Python><![CDATA[# Python code]]></Python>
    <Rust><![CDATA[# Rust code]]></Rust>
  </implementations>
  
  <propagation_status>
    <implementation lang="R" status="complete" date="2024-01-15"/>
    <implementation lang="Python" status="pending"/>
  </propagation_status>
</golden_test>
```

## Directory Structure

```
golden_tests/
├── specs/
│   ├── core/
│   │   ├── event_model/
│   │   │   ├── basic_hrf.xml          # Start here
│   │   │   ├── multiple_conditions.xml
│   │   │   └── continuous_regressors.xml
│   │   ├── baseline_model/
│   │   └── hrf_bases/
│   ├── integration/
│   └── edge_cases/
├── schema/
│   └── golden_test.xsd
└── validators/
    ├── R/validate_specs.R
    ├── Python/validate_specs.py
    └── Rust/validate_specs.rs
```

## Validation Types

1. **Dimensional**: Matrix/array dimensions, shape consistency
2. **Value Checks**:
   - Exact: For integers/categorical
   - Approximate: Floating-point with tolerance
   - Range: Within bounds
   - Relative: Percentage-based
3. **Statistical**: Sum, mean, std dev, min/max, percentiles
4. **Structural**: Column/row names, ordering, sparsity

## Workflow for New Tests

1. **Create XML spec** in appropriate directory
2. **Write semantic description** (purpose + algorithm)
3. **Implement in R** and validate outputs
4. **Document propagation status**
5. **Other languages** see failing tests and implement
6. **Update XML** with their implementations

## Workflow for New Language Implementation

1. **Get test specs** (submodule/copy from R repo)
2. **Create validator** in your test suite:
   ```python
   class GoldenTestValidator:
       def parse_golden_test(xml_path)
       def perform_numeric_check(matrix, check)
       def validate_test(xml_path)
   ```
3. **Start with basic_hrf.xml** - simplest test
4. **Implement required functionality** based on semantic descriptions
5. **Add your code to XML** implementations section
6. **Submit PR** to share implementation

## Test Sharing Methods

1. **Git Submodule** (recommended):
   ```bash
   git submodule add https://github.com/user/fmridesign.git fmridesign-r
   ln -s fmridesign-r/golden_tests golden_tests
   ```

2. **Separate Repository**: Dedicated golden-tests repo

3. **Package Distribution**: Include in package data

## Best Practices

### Adding Tests
- One behavior per test
- Start simple, minimal inputs
- Document WHY test exists
- Set appropriate tolerances (tighter for deterministic, looser for iterative)

### Implementing in New Language
- Match semantics, not syntax
- Focus on numeric equivalence
- Document any deviations in `<implementation_notes>`
- Use idiomatic code for your language

### Handling Differences
```xml
<implementation_notes>
  <note lang="Python">
    Uses scipy.linalg, tolerance 1e-6 for eigenvalues
  </note>
</implementation_notes>
```

## Current Test Status

| Test ID | R | Python | Rust |
|---------|---|--------|------|
| event_model_basic_hrf | ✅ | ⏳ | ⏳ |
| event_model_multiple_conditions | ✅ | ⏳ | ⏳ |
| baseline_model_polynomial_drift | ✅ | ⏳ | ⏳ |

Legend: ✅ Complete, ⏳ Pending, 🚧 In Progress

## Key Files to Read for Context

1. **Test examples**: Look at `specs/core/event_model/basic_hrf.xml`
2. **Schema**: Review `schema/golden_test.xsd` for XML structure
3. **Validators**: Check language-specific validators for implementation patterns

## Critical Points for AI Understanding

1. **Tests define behavior, not implementation** - focus on mathematical equivalence
2. **All validation is numeric** - no string comparisons, UI testing, or performance metrics
3. **Each language maintains own validator** - not in golden_tests directory
4. **Progressive workflow** - R implements first, others follow semantic spec
5. **Tolerances matter** - document why specific values chosen
6. **Cross-language collaboration** - implementations shared via XML updates

## Common Pitfalls to Avoid

- Don't modify existing specs when adding new language
- Don't test implementation details (data structures, variable names)
- Don't assume specific libraries available
- Don't use language-specific features in semantic descriptions
- Don't forget to update propagation_status

## Lessons Learned and Best Practices

### 1. Always Execute Code Before Writing Tests
- **Principle**: Never assume API behavior - verify through execution
- **Why it matters**: Function signatures, return types, and data structures often differ from expectations
- **Best practice**: 
  - Run code interactively before writing test specifications
  - Verify actual outputs match your mental model
  - Document any surprising behaviors in test comments

### 2. Be Aware of Function Polymorphism
- **Principle**: Many functions have multiple signatures with different behaviors
- **Why it matters**: The same function name may process arguments differently based on type or count
- **Best practice**:
  - Test all relevant function signatures
  - Read documentation for overloaded methods
  - Example: A function might treat `func(vec)` vs `func(x, y, z)` completely differently

### 3. Handle Object Systems Appropriately
- **Principle**: Different languages use different object models (S3/S4/R6 in R, classes in Python, structs in Rust)
- **Why it matters**: Direct field access, type coercion, and method calls vary by system
- **Best practice**:
  - Use appropriate accessor methods rather than direct field access
  - Test type conversions explicitly
  - Don't assume automatic coercion will work

### 4. XML Encoding Requirements
- **Principle**: XML has reserved characters that must be escaped
- **Common escapes**:
  - `<` → `&lt;`
  - `>` → `&gt;`
  - `&` → `&amp;`
  - `"` → `&quot;`
  - `'` → `&apos;`
- **Best practice**: Always escape comparison operators and special characters in test expressions

### 5. Verify Data Structure Internals
- **Principle**: Don't assume field names or structure without verification
- **Why it matters**: Internal representations often differ from external documentation
- **Best practice**:
  - Inspect objects programmatically (e.g., `str()` in R, `dir()` in Python)
  - Check field names exactly as they appear
  - Verify nested structure assumptions

### 6. Account for Language-Specific Conventions
- **Common differences**:
  - Indexing: 0-based (Python, Rust) vs 1-based (R, MATLAB)
  - Naming: snake_case vs camelCase conventions
  - Parameters: Positional vs named arguments
  - Types: Static vs dynamic typing implications
- **Best practice**: Document these differences in test comments when relevant

### 7. Recommended Test Development Workflow
1. **Prototype**: Write minimal working code in target language
2. **Execute**: Run code and capture actual outputs
3. **Inspect**: Examine data types, structures, and edge cases
4. **Specify**: Write golden test XML based on verified behavior
5. **Validate**: Run tests immediately to catch specification errors
6. **Iterate**: Refine based on validation results

### 8. Debugging Strategies for Failed Tests
- **Enable verbose output**: Most validators have detailed debugging modes
- **Isolate failures**: Extract failing checks into standalone scripts
- **Compare systematically**:
  - Expected vs actual values
  - Data types and structures
  - Numeric precision issues
- **Common issues**:
  - Floating-point comparison without tolerance
  - Wrong accessor methods for complex objects
  - Unescaped XML characters
  - Incorrect array/matrix dimensions

This consolidated reference provides everything needed to understand and work with golden tests efficiently.