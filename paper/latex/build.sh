#!/usr/bin/env bash
# Build the certified-atlas preprint from the markdown source, house template:
#   1. pandoc md -> body.tex (title block rendered by main.tex, so the body
#      starts at the italic preprint note)
#   2. python post-process: full-width rules, [N] -> \cite{N}, the markdown
#      References section replaced by \input of the 3.5-style bibliography,
#      the "Figures" list turned into real floats placed at first mention,
#      ToC insertion, breakable \texttt
#   3. pdflatex x2 for ToC and cross-references
# main.tex, preamble.tex and the bibliography are hand-maintained; the
# body is regenerated on every run, and the markdown is the SINGLE source of
# record for the text.
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"
PANDOC="${PANDOC:-$(command -v pandoc)}"
SRC=../certified_k3_atlas.md

python3 - <<'PY'
import re
lines = open('../certified_k3_atlas.md', encoding='utf-8').read().split('\n')
start = next(i for i, ln in enumerate(lines) if ln.startswith('*Preprint.'))
body = [re.sub(r'\[([^\]]+)\]\(#[^)]+\)', r'\1', ln) for ln in lines[start:]]
open('_body.md', 'w', encoding='utf-8').write('\n'.join(body))
PY

"$PANDOC" _body.md -f markdown -t latex --wrap=none --no-highlight \
    --top-level-division=section --shift-heading-level-by=-1 -o certified_k3_atlas_body.tex

python3 - <<'PY'
import re, os
p = 'certified_k3_atlas_body.tex'
txt = open(p, encoding='utf-8').read()

# --- full-width horizontal rules (3.5 template) ---
txt = txt.replace(r'\rule{0.5\linewidth}', r'\rule{\linewidth}')

# --- citations: [N], [N, M], [N-M] -> \cite{...} (numbered 3.5 bibliography) ---
NREF = 17
def cite(m):
    inner = m.group(1)
    out = []
    for part in inner.split(','):
        part = part.strip()
        rng = re.fullmatch(r'(\d+)\s*(?:-|--|–)\s*(\d+)', part)
        if rng:
            a, b = int(rng.group(1)), int(rng.group(2))
            if not (1 <= a <= b <= NREF):
                return m.group(0)
            out += [str(k) for k in range(a, b + 1)]
        elif part.isdigit() and 1 <= int(part) <= NREF:
            out.append(part)
        else:
            return m.group(0)
    return r'\cite{' + ','.join(out) + '}'
# pandoc escapes brackets as {[} ... {]} in text
txt = re.sub(r'\{\[\}([0-9,\s\-–]+)\{\]\}', cite, txt)
txt = re.sub(r'(?<!\\cite)\[([0-9,\s\-–]+)\]', cite, txt)

# --- the markdown References section becomes the 3.5 bibliography ---
m = re.search(r'\\section\{References\}(.*?)(?=\\section\{Appendix A)', txt, re.S)
assert m, 'References section not found; sync build.sh'
txt = txt[:m.start()] + '\\input{certified_k3_atlas_bib}\n\n' + txt[m.end():]

# --- figures: the markdown "Figures" list becomes real floats at first mention ---
fm = re.search(r'\\section\{Figures\}(.*)$', txt, re.S)
assert fm, 'Figures section not found; sync build.sh'
figblock = fm.group(1)
figs = re.findall(r'\\item\s*\n?\s*\\texttt\{([^}]*?\.png)\}\s*(.*?)(?=\n\s*\\item|\Z)', figblock, re.S)
if not figs:
    figs = re.findall(r'\\texttt\{([A-Za-z0-9_\\]*?\.png)\}\s*(.*?)(?=\\texttt\{[A-Za-z0-9_\\]*?\.png\}|\Z)', figblock, re.S)
assert len(figs) == 5, f'expected 5 figures, parsed {len(figs)}; sync build.sh'
txt = txt[:fm.start()]                      # drop the list section itself

def clean_caption(c):
    c = re.sub(r'\\(begin|end)\{[^}]*\}', ' ', c)
    c = c.replace(r'\item', ' ').replace('\n', ' ')
    c = re.sub(r'^\s*-{1,2}\s*', '', c.strip())
    c = re.sub(r'\s+', ' ', c).strip().rstrip('.')
    return c[:1].upper() + c[1:] if c else c

for n, (fn, cap) in enumerate(figs, start=1):
    name = fn.replace('\\_', '_')
    caption = clean_caption(cap)
    float_tex = ('\n\n\\begin{figure}[htbp]\n\\centering\n'
                 f'\\includegraphics[width=0.86\\linewidth,height=0.42\\textheight,keepaspectratio]{{{name}}}\n'
                 f'\\caption{{{caption}}}\n'
                 f'\\label{{fig:{n}}}\n\\end{{figure}}\n')
    # insert after the paragraph that first mentions "Figure n"
    hit = re.search(r'Figure[~\s]' + str(n) + r'\b', txt)   # NB: the space is required here ("Figure 1")
    if hit:
        nxt = txt.find('\n\n', hit.end())
        pos = nxt if nxt != -1 else len(txt)
        txt = txt[:pos] + float_tex + txt[pos:]
    else:
        txt += float_tex

# --- breakable long identifiers inside \texttt (3.5 template) ---
def breakable_tt(m):
    inner = m.group(1)
    inner = re.sub(r'(?<=[a-z])(?=[A-Z])', r'\\allowbreak{}', inner)
    inner = inner.replace(r'\_', r'\_\allowbreak{}')
    inner = re.sub(r'([./-])', r'\1\\allowbreak{}', inner)
    return r'\texttt{' + inner + '}'
txt = re.sub(r'\\texttt\{([^{}]*)\}', breakable_tt, txt)

# --- ToC between the abstract and the body proper ---
# The table of contents must NOT inherit the body's \parskip: at 0.4em per
# entry that costs about five lines over nineteen entries, which was enough to
# push Appendix F onto a page of its own.
TOC = ('\\newpage\n'
       '\\begingroup\n'
       '\\setlength{\\parskip}{0pt}\n'
       '\\setstretch{1.0}\n'
       '\\tableofcontents\n'
       '\\endgroup\n'
       '\\newpage\n')
txt = txt.replace(r'\section{1. Introduction}', TOC + '\\section{1. Introduction}', 1)

open(p, 'w', encoding='utf-8').write(txt)
print(f'  post-processed: {len(figs)} figures, citations converted')
PY

pdflatex -interaction=nonstopmode certified_k3_atlas.tex > _build.log 2>&1 || true
pdflatex -interaction=nonstopmode certified_k3_atlas.tex > _build.log 2>&1 || true
if [ -f certified_k3_atlas.pdf ]; then
  pp="$(pdfinfo certified_k3_atlas.pdf 2>/dev/null | awk '/Pages/{print $2}')"
  err="$(grep -ac '^!' _build.log || true)"
  und="$(grep -ac 'undefined' _build.log || true)"
  echo "  OK  certified_k3_atlas.pdf  (${pp} pp, ${err} errors, ${und} undefined-ref warnings)"
else
  echo "  FAILED (see _build.log)"; exit 1
fi
