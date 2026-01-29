# C13NV


```@eval
using Markdown
using Pkg

VERSION = Pkg.dependencies()[Base.UUID("136d2b59-c380-4582-9748-dc44afbd1b80")].version

github_badge = "[![Github](https://img.shields.io/badge/ARLQCI-C13NV.jl-blue.svg?logo=github)](https://github.com/ARLQCI/C13NV.jl)"

version_badge = "![v$VERSION](https://img.shields.io/badge/version-v$(replace("$VERSION", "-" => "--"))-green.svg)"

Markdown.parse("$github_badge $version_badge")
```

## Contents

```@contents
Depth = 2
Pages = [pair[2] for pair in Main.PAGES[2:end-1]]
```
