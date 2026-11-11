# Evidence record: Solyanik estimates for the strong maximal Tauberian function

**Checked:** 2026-11-11
**Question:** Does the published Solyanik estimate, with its actual dimension constants and near-one range, imply or contradict a bound along \(n=N\to\infty\), \(\alpha_N=1-\tau_N-B\log N/N\), with \(\tau_N\sim\sqrt{\mu\log N/N}\)?

## Sources and scope

### S1 — A. A. Solyanik, original theorem

- **Identity:** A. A. Solyanik, “On halo functions for differentiation bases,” *Matematicheskie Zametki* 54(6) (1993), 82–89; English translation, *Mathematical Notes* 54(5–6) (1993/1994), 1241–1245. DOI `10.1007/BF01209086`; MR1268374. MathNet's English catalog renders the Russian title as “Bordering functions for differential bases.”
- **Primary record/full text:** [MathNet record](https://www.mathnet.ru/eng/mzm2453), [MathNet full-text PDF](https://www.mathnet.ru/php/getFT.phtml?jrnid=mzm&option_lang=eng&paperid=2453&what=fullt).
- **Type/version:** Original journal article, read in the original Russian scan; bibliographic record identifies the English translation.
- **Read scope:** PDF pages 2–6 (printed pages 83–87), mainly printed p.83 theorem/notation, pp.85–87 proof; record metadata also checked. This is a localized reading, not a line-by-line transcription of all eight pages.
- **Relevant statement:** For the strong maximal function over axis-parallel rectangles, the original defines its halo function \(\phi_n(u)\). For \(n>2\), it states \(\phi_n(u)-1\asymp (u-1)^{1/n}\) as \(u\downarrow1\). The paper defines its near-one comparison notation by positive constants bounding the ratio on \(1<u<2\). Thus the theorem supplies two-sided bounds on that range, but the constants are not numerically specified in the statement.
- **Excerpt (10 words):** “for n>2, φ_n(u)−1 ∼ (u−1)^{1/n}” (formula counted conservatively).
- **Limit:** The notation on the theorem is dimension indexed; no dimension-uniform values for its lower/upper constants are given. The proof introduces dimension-dependent induction constants and does not provide estimates usable here as an absolute constant independent of \(n\).

### S2 — Hagelstein–Parissis, English restatement and dimensional notation

- **Identity:** Paul Hagelstein and Ioannis Parissis, “Solyanik estimates in ergodic theory,” *Colloquium Mathematicum* 145(2) (2016), 193–207. DOI `10.4064/cm6729-2-2016`.
- **Publisher record/full text:** [IMPAN article record](https://www.impan.pl/en/publishing-house/journals-and-series/colloquium-mathematicum/all/145/2/91594/solyanik-estimates-in-ergodic-theory), [publisher PDF](https://www.impan.pl/shop/en/publication/transaction/download/product/91594).
- **Type/version:** Peer-reviewed journal article, publisher-hosted final PDF; used here as a secondary statement of the geometric theorem (the paper's main new results are ergodic analogues).
- **Read scope:** PDF pages 2–4 (printed pages 194–196), definitions and Theorem 1.3; PDF page 13 (printed p.205), notation; references on PDF pages 14–15.
- **Relevant statement:** The paper defines \(C_S(\alpha)\) as the sharp restricted Tauberian constant for the uncentered strong maximal operator over axis-parallel rectangles in \(\mathbb R^n\). Theorem 1.3 states \(C_S(\alpha)-1\lesssim_n (1/\alpha-1)^{1/n}\) for \(\alpha\) sufficiently close to 1, and cites Hagelstein–Parissis and Solyanik. The notation paragraph explicitly says the subscript \(n\) denotes dependence of the implicit constant on dimension.
- **Excerpt (7 words):** “for α sufficiently close to 1”.
- **Limit:** This is not a new uniform-in-dimension theorem. It gives no explicit numerical value for either the multiplicative constant or the “sufficiently close” threshold as a function of \(n\).

### S3 — Hagelstein–Parissis, harmonic-analysis exposition

- **Identity:** Paul Hagelstein and Ioannis Parissis, “Solyanik estimates in harmonic analysis,” in *Special Functions, Partial Differential Equations, and Harmonic Analysis*, Springer Proceedings in Mathematics & Statistics 108 (Springer, 2014), 87–103. DOI `10.1007/978-3-319-10545-1_9`.
- **Full text:** [author manuscript/preprint](https://arxiv.org/abs/1310.3771), [opened PDF copy](https://citeseerx.ist.psu.edu/document?doi=fae3e600f08aa9ca11d3d5b1229d27d104e588a1&repid=rep1&type=pdf).
- **Type/version:** Published proceedings chapter; read the available full-text manuscript/PDF. This is a secondary restatement with proof discussion, not the original theorem.
- **Read scope:** PDF/printed pp.4–9 for Tauberian definition, Solyanik theorem, notation and the explicit estimate for the iterated coordinate maximal operator; references pp.17–18.
- **Relevant statement:** It records \(C_S(\alpha)-1\sim_n (1/\alpha-1)^{1/n}\) as \(\alpha\to1^-\), attributes the strong estimate to Solyanik, and defines \(\sim_n\) through two inequalities whose implicit constants may depend on \(n\).
- **Excerpt (9 words):** “the constant c depends … on the dimension n”.
- **Additional explicit bound in this source:** Printed Theorem 2 gives, for the iterated one-coordinate maximal operator \(M_1\cdots M_n\),
  \[
  C_{1\cdots n}(\alpha)\le
  \left(1+4\frac{(1-\alpha)^{1/n}}{1-(1-\alpha)^{1/n}}\right)^n,
  \qquad 0<\alpha<1.
  \]
  Since each axis-parallel rectangle average is bounded by the corresponding iterated coordinate maximal expression, this is also a valid explicit (coarse) upper bound for the geometric strong Tauberian constant \(C_S(\alpha)\). This is not the implicit constant in the sharp \(\sim_n\) asymptotic, and should not be described as an absolute-constant Solyanik estimate.
- **Limit:** The sharp asymptotic comparison itself has dimension-dependent constants. The explicit all-\(\alpha\) upper bound above is available, but along \(\varepsilon_N=1-\alpha_N\asymp\sqrt{\log N/N}\) its logarithm is of order at most \(N\log(N/\log N)\) by direct substitution, much weaker than the proposed \(\mu N+O(\sqrt{N\log N})\) scale. It does not settle that target or contradict an exponential lower bound of that scale.

## Interpretation for the proposed sequence

Set \(u=1/\alpha\). Then the original theorem has \(u\in(1,2)\), corresponding to \(\alpha\in(1/2,1)\). The sequence in question has

\[
\varepsilon_N=1-\alpha_N=\tau_N+B\frac{\log N}{N}\to0,
\]

so for every fixed finite \(B\), it is eventually in this near-one range. In terms of \(\varepsilon_N\), the theorem gives, for each fixed dimension \(N>2\),

\[
c_N\left(\frac{\varepsilon_N}{1-\varepsilon_N}\right)^{1/N}
\le H_N(\alpha_N)-1\le
C_N\left(\frac{\varepsilon_N}{1-\varepsilon_N}\right)^{1/N},
\]

with positive constants \(c_N,C_N\) depending on dimension, and with no usable dimension-uniform values supplied by these sources. The power term is

\[
\left(\frac{\varepsilon_N}{1-\varepsilon_N}\right)^{1/N}
=\exp\!\left(\frac{\log\varepsilon_N-\log(1-\varepsilon_N)}{N}\right)\longrightarrow1.
\]

This confirms applicability of the *fixed-dimension theorem's stated \(\alpha\)-range* for sufficiently large \(N\), but the sharp asymptotic constants \(C_N,c_N\) are dimension dependent. The published explicit all-\(\alpha\) bound for an operator dominating the geometric strong maximal operator is available, but at this sequence it is much coarser, with logarithm on the \(N\log(N/\log N)\) scale. Neither the sharp asymptotic statement nor this explicit coarse estimate yields the requested \(\log H_N\le\mu N-d_N+c_Q\log N+o(\log N)\). Neither supplies a contrary bound strong enough to conflict with an internal exponential lower bound of scale \(\mu N-d_N\). The literature reviewed here therefore reveals no contradiction. No numerical value such as \(c_Q\approx0.3217\) appears in these Solyanik results.

## Verification and bibliographic checks

- S1 author/title/venue/year/volume/issue/pages/translation pages/DOI were checked against MathNet's journal record and its AMSBIB entry; the author-hosted/later bibliographies agree.
- S2 author/title/venue/year/volume/issue/pages/DOI were checked on the publisher's IMPAN article record; theorem and notation were read in its publisher PDF.
- S3 title/authors/chapter pages/series/volume/publisher/year/DOI were checked on the Springer book page and the chapter manuscript/PDF metadata; the result is included as a useful secondary cross-check, not represented as an original theorem.
- Search date and actual query strings are recorded in `search-log.md`. Bibliographic records are limited to these three sources; no unverified entries are included.

## Parent verification and metadata correction

The parent opened the MathNet record, the IMPAN publisher record, arXiv v2 and the Springer book/chapter pages on 2026-11-11. The published book subtitle is **In Honor of Calixto P. Calderón**, corrected in the BibTeX; its chapter DOI is the single-underscore `10.1007/978-3-319-10545-1_9`, confirmed by the chapter’s publisher citation. The arXiv metadata displays a double underscore, so that rendering was not exported. The Russian original is 1993. MathNet catalogs the translation as 54(6), 1993; Springer’s reference gives 54(5–6), 1994. This translation catalog discrepancy does not change the original publication date and is retained here rather than silently merging years.

Primary checks: [Springer book](https://link.springer.com/book/10.1007/978-3-319-10545-1), [Springer chapter](https://link.springer.com/chapter/10.1007/978-3-319-10545-1_9), [author manuscript PDF](https://arxiv.org/pdf/1310.3771). No priority claim for the present internal results is made.
