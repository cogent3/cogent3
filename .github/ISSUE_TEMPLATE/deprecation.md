---
name: Flag for deprecation
about: identify a function / method / attribute for deprecation
title: ''
labels: ''
assignees: ''

---

**What should be changed**

Deprecated means a function, method or argument is going to be replaced in the code base.

What function, method or argument should be marked for change?

What is the full module path to the attribute?

**Why**

What should it be replaced with and why?

**How much notice**

**How much notice**

How much time should the community be given to response to this notice? Specify as the number of months (*x*) before the capability will be removed. Leave blank if not known.

**Guidelines for deprecation**

As this is a function that will be removed, see [dev guidelines for deprecation](https://github.com/cogent3/c3dev/wiki/Deprecating-Code) of code. Indicate the limiting version (for when the old function will be removed) as <current year>.<PR month + *x*> where *x* is defined above.
