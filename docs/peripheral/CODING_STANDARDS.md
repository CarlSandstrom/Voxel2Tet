# ESD Tool Development Code Standards (BEAST) (Orpheus), version 1.1

## Table of Contents
- Review of code standards
- Introduction
- Scope
- Portability
- Naming
- Code Design
- Formatting
- Documentation

---

## Review of code standards

## Introduction

This document describes code standards that apply to all C++ development work within ESD. The code standards are based on:

- Ellemtel-standards
- *Effective C++* by Scott Meyers
- *C++ Coding Standards*
- Common practice within the BEAST and Orpheus teams

## Scope

The rules apply to all C++ code developed within ESD. Exceptions may apply for specific tools and are explicitly stated.

## Portability

Compilers used within ESD support **C++17**, but not higher. Do not use features from newer standards.

---

## Naming

### Styles

- `CamelCase` – Capitalize each word, no underscores or digits.
- `camelCase` – Like CamelCase but starts with lowercase.
- `ALL_CAPITALS` – All caps, words separated by underscores.

### General

- Use English for all identifiers.
- Use clear, unabbreviated names unless abbreviation is well-established (e.g., CAD).
- Names should reflect purpose.

### Directories

- Use descriptive `CamelCase` names.

### Files

- Match file name to class name (`CamelCase`).
- Use `.h` for headers, `.cpp` for implementations  
  _BEAST Exception: Use `.cc` instead._
- Use `.tpp` for template implementations.

### Classes

- Use `CamelCase`.
- Name should reflect purpose.

### Functions

- Use `camelCase`.
- Setters start with `set`, getters with `get`.
- Mathematical notation may use underscores (e.g., `get_dp_dx()`).

### Variables

- Non-class: `camelCase` (e.g., `wearRate`)
- Class members: `m_` prefix (e.g., `m_wearRate`)  
  _BEAST Exception: Use trailing underscore (e.g., `wearRate_`)_
- Static members: `ms_` prefix  
  _BEAST Exception: Ignore this rule._
- Namespace variables: `n_` prefix  
  _BEAST Exception: Ignore this rule._
- Avoid magic numbers.
- Constants: `ALL_CAPITALS`
- Boolean: Start with `is`, `has`, `should` (e.g., `m_isInitialized`)
- Object instances should reflect class name:
  ```cpp
  Animal* animal;
  Dog niceDog;
  Dog badDog;
  ```

### Enums

- Enum class names: `CamelCase`
- Enum values: `ALL_CAPITALS`

### Structures

- Struct names: `CamelCase`
- Members follow class member naming conventions

### Namespaces

- Use `CamelCase`

### Typedefs

- Format: `CamelCase_t`
  ```cpp
  typedef std::map<std::string, Lubricant> LubricantMap_t;
  ```

---

## Code Design

### General

- Only use `#pragma once`  
  _BEAST Exception: `#ifdef` allowed for debug toggles._
- Resolve compiler warnings before pushing code.

**GCC Flags:**
```text
-O3
-Weffc++
-Wall
-Wextra
-Wshadow
-Werror
-Wc++0x-compat
-Wcast-qual
-Wlogical-op
-Wnon-virtual-dtor
-Wold-style-cast
-Woverloaded-virtual
-fdiagnostics-show-option
```

### Classes

#### General

- One class per file.
- No public member variables.
- Friend classes only for unit testing.
- Follow the "rule of five":
  - Destructor
  - Copy constructor
  - Copy assignment
  - Move constructor
  - Move assignment
- Use `final` to prevent inheritance.
- Use `delete` to prevent instantiation.
- Use `default` for default implementations.
- Do not change access via inheritance.

#### Constructors & Destructors

- Base classes must have virtual destructors.
- Use initializer lists.
- Use `explicit` for single-argument constructors  
  _BEAST Exception: Ignore this rule._

#### Assignment Operator

Example:
```cpp
B& B::operator=(const B& source) {
  if (this == &source) return *this;
  A::operator=(source);
  // assignments
  return *this;
}
```

- Prevent self-assignment.
- Return reference.
- Assign all members, including base class.

### Functions

- Pass non-primitive types by `const` reference or pointer.
- Use `const` for non-modifying member functions.
- Do not return non-`const` references from `const` functions.
- Avoid non-`const` pointers/references unless documented.
- Avoid non-`const` member pointers.
- Implement functions consistently (all in header or implementation).
- Use `override` and `final` appropriately.
- Do not redefine inherited non-virtual functions or default parameters.
- Avoid calling virtual functions in constructors/destructors.
- Maintain access specifiers across inheritance.
- Getters: `const`, Setters: `void`
- Getters/setters should be simple.

### Variables

- No public/protected member data.
- Use `const` and `constexpr` appropriately.
- Prefer `enum class` over `enum`.
- Declare variables in smallest scope.
- Avoid `mutable`.
- Use `nullptr`.
- Prefer pointers over references.
- Use `unique_ptr`, avoid `shared_ptr` unless necessary.
- Never use `auto_ptr`.
- Define static members in `.cpp`.
- Avoid global variables.

### Exceptions

- Choose one documented exception protocol:
  1. Document exceptions per function (`@throws`)
  2. Document globally (e.g., in README)
- Use custom exception base class derived from `std::exception`.
- Throw early on invalid state.
- Avoid throwing in destructors.
- Catch exceptions by `const &`.
- Avoid `noexcept`.

### Statements

- Only use `static_cast` and `dynamic_cast`.
- Cast to highest necessary type.
- `switch` must have `default`.
- Avoid `goto`.
- Avoid manual memory management (`new`, `delete`, etc.).
- Prefer `break`, `return`, `continue` over deep nesting.
- Prefer range-based loops and STL algorithms.
- Use `assert()` for development errors.
- Use exceptions for runtime errors.
- Use `auto` where appropriate.

### Namespaces

- Do not use `using namespace` in headers.
- In `.cpp`, either use `using namespace`, wrap code, or use specific `using`.

### Header Files

- Use `#pragma once`.
- Prefer forward declarations over includes.
- Never include files outside project directory.
- Use `""` for project includes, `< >` for Qt/standard.
- Include order:
  1. Base class
  2. Project headers
  3. Third-party
  4. Standard library
- Forward declarations follow includes.
- Template `.tpp` files go at bottom of header.

### Implementation Files

- Include only relevant headers.
- Include order:
  1. Corresponding header
  2. Project headers
  3. Third-party
  4. Standard library

---

## Formatting

### Manual Formatting

```cpp
void myFunction() {
  ...
  return;
}
```

- One variable per line.
- Class member order:
  1. Enums/constants
  2. Functions
  3. Variables (only private)
- Qt signals/slots order:
  - public
  - public slots
  - signals
  - protected
  - protected slots
  - private (functions)
  - private slots
  - private (variables)
- Max 120 characters per line.
- Max 5 indentation levels.
- Max 20 lines per function.
- No spaces around `.` or `->`.
- No space between function name and brackets.

### Automatic Formatting (Beautifier)

- Uses `astyle` or SKF rules  
  _BEAST Exception: Uses Clang formatter._
- Examples of rules:
  - `--style=bsd`: Braces on new lines
  - `--indent=spaces=2`  
    _BEAST Exception: Use 4 spaces_
  - `--break-blocks`: Pad with empty lines
  - `--add-braces`: Add braces to one-line conditionals  
    _BEAST Exception: Deactivated_
  - `--pad-paren-in`: Space inside parentheses  
    _BEAST Exception: Deactivated_
  - `--pad-oper`: Space around operators
  - `--convert-tabs`: Convert tabs to spaces
  - `--align-pointer=type`: Align `*` with type
  - `--align-reference=type`: Align `&` with type
  - `EmptyLineFormatter`: Remove duplicate empty lines
  - `IncludeGuardFormatter`: Replace include guards with `#pragma once`
  - `ExtraSpaceAroundBracketsFormatter`: Remove extra space in templates

---

## Documentation

### General

- All comments in English.
- Use **Doxygen** for documentation.
- Use `/** */` for blocks:
  ```cpp
  /**
   * @brief sum
   * @param a the first term
   * @param b the second term
   * @return the sum
   */
  int sum(int a, int b) {
    return a + b;
  }
  ```
- Use `///<` for member variable comments:
  ```cpp
  private:
    m_myMemberVariable ///< description of the variable
  ```

---
