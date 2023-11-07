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

The number of months (*x*) long before the capability will be removed? This should be represented in the user warning as a specific version. e.g. <current year>.<current month + *x*> where *x* is from the PR date.  

**Guidelines for deprecation**

As this is a function being deprecated, see [dev guidelines for deprecation](https://github.com/cogent3/c3dev/wiki/Deprecating-Code) of code. Indicate the limiting version (for when the old function will be removed) as 3 months from when you make the change.

