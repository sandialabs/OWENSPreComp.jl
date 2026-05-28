```@meta
CurrentModule = OWENSPreComp
```

# API Reference

The exported API is intentionally small. The tested workflow also relies on
module-qualified constructors and legacy reader helpers, so they are included
here for discoverability.

## Contents

```@contents
Pages = ["reference.md"]
Depth = 3
```

## Public Workflow

The main user-facing objects are `properties`, `tw_rate`,
`OWENSPreComp.Input`, and `OWENSPreComp.Output`. `properties` and `tw_rate` are
exported; the types are currently accessed with module qualification.

## Kernel Types And Functions

```@autodocs
Modules = [OWENSPreComp]
Pages = ["main.jl"]
Order = [:type, :function]
```

## Legacy File Readers

```@autodocs
Modules = [OWENSPreComp]
Pages = ["io.jl"]
Order = [:function]
```

## Deprecated Compatibility Surface

```@autodocs
Modules = [OWENSPreComp]
Pages = ["deprecated.jl"]
Order = [:type, :function]
```

## Complete Index

```@index
```
