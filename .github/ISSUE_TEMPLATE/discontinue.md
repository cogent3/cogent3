---
name: Flag for discontinued
about: Suggest an idea for this project
title: ''
labels: ''
assignees: ''

---

**What should be changed**

Discontinued means a function, method or argument is going to be deleted entirely from the code base and is not being replaced. 

What function, method or argument should be marked for removal?

What is the full module path to the attribute?

**Why**

Why should it be discontinued?

**How much notice**

The number of months (*x*) long before the capability will be removed? This should be represented in the user warning as a specific version. e.g. <current year>.<current month + *x*> where *x* is from the PR date.  

**Guidelines for deprecation**

As this is a function that will be removed, see [dev guidelines for deprecation](https://github.com/cogent3/c3dev/wiki/Deprecating-Code) of code. Indicate the limiting version (for when the old function will be removed) as 3 months from when you make the change.

