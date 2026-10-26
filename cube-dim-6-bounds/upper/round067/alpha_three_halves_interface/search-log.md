# 067 primary-source reading log

Date: 2026-10-26. Narrow goal: verify the Fourier normalization, scaling and one-dimensional stable-density bound used by the already existing alpha=3/2 body32 interface. This was not an exhaustive literature review. No new stable-density computation was run.

Actual web search queries:

1. `symmetric alpha stable density estimate t |x| 1 alpha Chen Kim Song heat kernel estimates stable processes arxiv`
2. `stable densities symmetric alpha stable scaling tail p_t(x) t t^{1/alpha} primary paper`

The attempted direct DOI `https://doi.org/10.4171/JEMS/231` redirected to an unavailable EMS page. It is not counted as an independently read source and is not cited as support.

The successful source was Bogdan, Grzywny and Ryznar, *Heat kernel estimates for the fractional Laplacian with Dirichlet conditions*: [arXiv abstract/version page](https://arxiv.org/abs/0905.2626), [versioned PDF](https://arxiv.org/pdf/0905.2626v4), [published DOI](https://doi.org/10.1214/10-AOP532). Title, authors, journal reference, DOI and revision were checked against the abstract/PDF. Version v4 was revised 2010-11-05; the original submission was 2009-05-15.

Actual reading: metadata and introductory context; electronic-reprint PDF page 3, equation (4), and PDF page 5, equations (5)--(9). These support the isotropic stable Fourier normalization, probability convolution semigroup, scaling and two-sided density order. Only the d=1 specialization is used as corroboration for the independent coordinate factors. This does not transfer an isotropic six-dimensional estimate to body32. The whole paper and its earlier references were not independently reread.

The PDF was downloaded from the versioned URL using Python urllib. It contains 316572 bytes; SHA-256 `fe711b64b9de875b38aba1c95f25b0ce0c435eab981f1d671c6c72db1ae3e15b`. The actual download completed in session 30240, tool chunk 016f7e. It is retained under `sources/` because this task specifically requested version/hash binding.

Historical C006 sources.json contains Watanabe and Simon source records. Those records were read and hashed, but the cited papers were not independently reread in 067, and their historical download hashes are not labeled new 067 acquisitions. Five-dimensional transfer notes were read only as abstract background; no five-dimensional computation or quantitative status was inherited.

After discovery of C006/F006 and C007/F008/A008, the parent directed this branch to actual 0 science. The final work therefore audits provenance and accurately identifies the unpaid interior/competitive-pair problem; it does not repeat the old theoretical or numerical batch.
