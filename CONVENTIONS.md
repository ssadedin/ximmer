# Development Guide for Ximmer

- always use spaces for indentation and code formatting. Do NOT use tabs anywhere.
- assume a tab equivalent size of 4 characters for indentation
- for efficiency and correctness, strongly prefer static compilation in any place where
  truly dynamic behavior is not needed

NOTE: Although this project uses Groovy 3.x, it maintains 2.x syntax compatibility. Therefore, do not
use Groovy 3.x (Parrot Parser) syntax. For example:

- WRONG:  `new byte[] { 1,2,3 }`  : fails to compile on Groovy 2.x
- CORRECT: `[1,2,3] as byte[]`

## Quirks

- when accessing CLI options on CliOptions or OptionAccessor object, use the array style access so that
  static compilation works:
  - WRONG: `opts.foo`
  - CORRECT: `opts['foo']`

## Testing

- tests are implemented using junit4 exclusively
- for advanced mocking, use the Mockito library

It is PREFERRED to use native Groovy assertions in tests, rather than Junit assertions,
because they produce more expressive output.

Example:

```groovy
// GOOD - use native groovy assertion
assert foo.bar != 'baz'

// BAD - using Junit assertion when Groovy is more expressive
assertNotEquals(foo.bar, 'baz') // BAD
```
