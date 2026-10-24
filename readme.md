# Crosstalk Fortran Programs — Introduction to Electromagnetic Compatibility (Paul, Scully, Steffka, 3rd ed.)

Seven Fortran 77 programs by Clayton R. Paul (University of Kentucky) for computing crosstalk on Multiconductor Transmission Lines (MTLs). This README covers building them with **VS Code + WSL2** and explains what each source file does. The codes were take from the [Wiley FTP](https://ftp.wiley.com/public/sci_tech_med/electromagnetic_compatability) repository.

---

## Contents

1. [Package overview](#1-package-overview)
2. [Building with VS Code + WSL2](#2-building-with-vs-code--wsl2)
3. [Conventions common to all programs](#3-conventions-common-to-all-programs)
4. [Stage 1 — Per-unit-length L and C](#4-stage-1--per-unit-length-l-and-c)
5. [Stage 2 — SPICE models](#5-stage-2--spice-models)
6. [Known issues and caveats](#6-known-issues-and-caveats)
7. [Publishing to GitHub](#7-publishing-to-github)
8. [Copyright notice](#8-copyright-notice)

---

## 1. Package overview

The workflow has two stages:

```
Stage 1 — per-unit-length parameters (write PUL.DAT)
  WIDESEP   wires in a homogeneous medium (wide-separation approximation)
  RIBBON    ribbon cable (wires with dielectric insulation)
  PCB       lands on a board, no ground plane
  MSTRP     coupled microstrip (lands on a substrate over a ground plane)
  STRPLINE  coupled stripline (lands midway between two ground planes)

Stage 2 — SPICE subcircuit of the line (read PUL.DAT)
  SPICEMTL  exact lossless-line model (modal decoupling)
  SPICELPI  approximate single lumped-Pi model
```

Pick the cross-section, compute **L** and **C**, then generate a subcircuit to simulate crosstalk (NEXT/FEXT) in PSpice or LTspice.

| Source | Input | Output |
|---|---|---|
| `WIDESEP.FOR`  | `WIDESEP.IN`  | `PUL.DAT` |
| `PCB.FOR`      | `PCB.IN`      | `PUL.DAT` |
| `MSTRP.FOR`    | `MSTRP.IN`    | `PUL.DAT` |
| `RIBBON.FOR`   | `RIBBON.IN`   | `PUL.DAT` |
| `STRPLINE.FOR` | `STRPLINE.IN` | `PUL.DAT` |
| `SPICEMTL.FOR` | `SPICEMTL.IN` + `PUL.DAT` | `SPICEMTL.OUT` |
| `SPICELPI.FOR` | `SPICELPI.IN` + `PUL.DAT` | `SPICELPI.OUT` |

The shipped `XXX.OUT` files of the five Stage 1 programs are only *reference copies* of `PUL.DAT`. Those programs never write a `.OUT` file. The shipped `.EXE` files are old Windows/DOS binaries and do not run on Linux; build from the `.FOR` sources instead.

---

## 2. Building with VS Code + WSL2

### 2.1 Install the compiler in WSL2

In the Ubuntu (WSL) terminal:

```bash
sudo apt update
sudo apt install -y gfortran make
gfortran --version
```

### 2.2 Copy the files into the Linux file system

Work under `~` rather than `/mnt/c/...`. It is much faster and avoids permission issues.

```bash
mkdir -p ~/emc/cap9
cp "/mnt/c/Users/<your_user>/<path>/Crosstalk Fortran Programs/"* ~/emc/cap9/
cd ~/emc/cap9
```

Copy the `.IN`, `.OUT` and `PUL.DAT` files too. **Keep file names in UPPERCASE.** The code opens `'WIDESEP.IN'`, `'PUL.DAT'`, etc., and Linux file names are case-sensitive.

### 2.3 Open the folder in VS Code

1. On Windows, install the **WSL** extension (Microsoft) in VS Code.
2. From the WSL terminal, inside the folder, run `code .`. VS Code opens connected to WSL; "WSL: Ubuntu" appears in the bottom-left corner.
3. Optional: install the **Modern Fortran** extension *on the WSL side* for syntax highlighting. Its linter uses the installed `gfortran`.

### 2.4 Compile

Single program:

```bash
gfortran -std=legacy -O2 -o widesep WIDESEP.FOR
```

- The `.FOR` extension makes `gfortran` treat the file as fixed-form source (columns 7–72). Windows CRLF line endings are not a problem.
- `-std=legacy` suppresses warnings about the old features these programs use: `PAUSE 'Singular matrix'` in `GAUSSJ` and arithmetic IFs in `JACOBI`. Without it the code still compiles, with warnings.
- For debugging, use `-g -fcheck=all -Wall` instead of `-O2`.
- For extra precision (everything is single-precision `REAL`), add `-fdefault-real-8`.

To build all seven, create a `Makefile` in the folder:

```make
FC     = gfortran
FFLAGS = -std=legacy -O2
PROGS  = widesep pcb mstrp ribbon strpline spicemtl spicelpi

all: $(PROGS)

widesep:  WIDESEP.FOR  ; $(FC) $(FFLAGS) -o $@ $<
pcb:      PCB.FOR      ; $(FC) $(FFLAGS) -o $@ $<
mstrp:    MSTRP.FOR    ; $(FC) $(FFLAGS) -o $@ $<
ribbon:   RIBBON.FOR   ; $(FC) $(FFLAGS) -o $@ $<
strpline: STRPLINE.FOR ; $(FC) $(FFLAGS) -o $@ $<
spicemtl: SPICEMTL.FOR ; $(FC) $(FFLAGS) -o $@ $<
spicelpi: SPICELPI.FOR ; $(FC) $(FFLAGS) -o $@ $<

clean: ; rm -f $(PROGS)
```

Then run `make` in the VS Code integrated terminal (Ctrl+`).

To build with **Ctrl+Shift+B**, create `.vscode/tasks.json`:

```json
{
  "version": "2.0.0",
  "tasks": [{
    "label": "make",
    "type": "shell",
    "command": "make",
    "group": { "kind": "build", "isDefault": true },
    "problemMatcher": "$gcc"
  }]
}
```

### 2.5 Run and verify

The programs are not interactive. They read the `.IN` file and write to the current directory:

```bash
mkdir -p ref && cp *.OUT ref/                # back up the reference outputs first
./widesep                                    # reads WIDESEP.IN  -> writes PUL.DAT
diff --strip-trailing-cr PUL.DAT ref/WIDESEP.OUT   # should match the reference
./spicemtl                                   # reads SPICEMTL.IN + PUL.DAT -> SPICEMTL.OUT
```

Notes:

- **Run order:** `SPICEMTL` and `SPICELPI` read the `PUL.DAT` produced by one of the five Stage 1 programs. Run the cross-section program first.
- **Overwriting:** every run overwrites `PUL.DAT`. `SPICEMTL` and `SPICELPI` also **overwrite the reference `.OUT` files**, hence the backup above.
- **Array size (`MSIZE`):** `WIDESEP.FOR`, `SPICEMTL.FOR` and `SPICELPI.FOR` use `PARAMETER(MSIZE=2)`, i.e. at most 2 conductors plus the reference. Beyond that, `WIDESEP` writes only `INPUT DATA ERROR` and `SPICEMTL` stops without producing a model. Increase `MSIZE` in the source and recompile. In `PCB`, `MSTRP`, `STRPLINE` and `RIBBON` the limit is `N*NCDIV ≤ 100` (or `N*NF ≤ 100`).
- **`Singular matrix` message:** the program is paused inside `PAUSE`. Press Enter to continue, but the result is most likely wrong because of the input data.

---

## 3. Conventions common to all programs

**I/O.** Each program reopens units 5 and 6 (normally keyboard and screen) as files: `OPEN(UNIT=5,FILE='XXX.IN')` and `OPEN(UNIT=6,FILE='PUL.DAT')`. The `.IN` files are read with list-directed `READ(5,*)`, one value per line, in the fixed order given for each program below.

**Units.** `CMTM=2.54E-5` converts mils to metres. Almost all input dimensions are in **mils**. The exception is WIDESEP, which uses coordinates in metres for reference types 1 and 2.

**Constants.** `V0=2.997925E8` is the speed of light and `EPS = 1/(v0²·μ0)` is ε0. The programs rely on the homogeneous-medium relation

$$\mathbf{L} = \mu_0\varepsilon_0\,\mathbf{C}_0^{-1} = \mathbf{C}_0^{-1}/v_0^2$$

where **C₀** is the capacitance computed with the dielectric removed (εr = 1). Inductance does not depend on the dielectric, so every inhomogeneous-medium program solves the problem twice: without and with the dielectric.

**`PUL.DAT` format.** Upper triangle of the matrices, one element per line:

```
  1    1   3.12345E-07        =L(  1,  1)
  1    2   ...                =L(  1,  2)
  ...                         =C(...)   =C0(...)   (some programs)
NUMBER OF LANDS= ...          (echo of the input data at the end)
```

The SPICE programs read only `I, J, value` with list-directed input, ignore the `=L(…)` text and the rest of the file, and expect the L block first, then the C block. Do not change that order.

**Subroutine `GAUSSJ`** (Numerical Recipes). Gauss–Jordan elimination with full pivoting; returns **A⁻¹ in place in A**. The `IDENT` argument is the right-hand side and is also overwritten. This is harmless because the result in A does not depend on it, but after the first call `IDENT` is no longer the identity. `PAUSE 'Singular matrix'` means a singular matrix, almost always caused by invalid input geometry.

**Precision.** Everything is single-precision `REAL` (~7 digits), which is enough for the examples. With many subdivisions, consider `-fdefault-real-8`.

**Reference conductor.** For N+1 conductors the line matrices are N×N. In WIDESEP, MSTRP and STRPLINE the reference is a separate object (wire, plane or shield). In PCB and RIBBON the reference is **one of the lands/wires** (`IREF`), and the program reduces the N×N matrix to (N−1)×(N−1).

---

## 4. Stage 1 — Per-unit-length L and C

### 4.1 `WIDESEP.FOR` — wires, wide-separation approximation

Closed-form formulas, no numerical method. Valid when the spacing between wires is large compared with their radii, so the charge on each wire can be treated as uniform. Homogeneous medium (εr, μr).

**Input `WIDESEP.IN`:** `N`, `REFTYPE` (1 = wire, 2 = ground plane, 3 = cylindrical shield), `ER`, `MUR`, `R0` (radius of the reference wire or shield, mils), then for each wire: `NUM`, `R` (mils), `X`, `Y`.

The meaning of X and Y depends on `REFTYPE`:

| REFTYPE | X | Y | Implemented formulas |
|---|---|---|---|
| 1, wire | position (m) relative to wire 0 | position (m) | $l_{ii}=\frac{\mu}{2\pi}\ln\frac{d_{i0}^2}{r_{wi}r_{w0}}$, $l_{ij}=\frac{\mu}{2\pi}\ln\frac{d_{i0}d_{j0}}{d_{ij}r_{w0}}$ |
| 2, plane | height h (m) | horizontal position (m) | $l_{ii}=\frac{\mu}{2\pi}\ln\frac{2h_i}{r_{wi}}$, $l_{ij}=\frac{\mu}{4\pi}\ln\left(1+\frac{4h_ih_j}{s_{ij}^2}\right)$ |
| 3, shield | radial distance (mils) | angle (degrees) | $l_{ii}=\frac{\mu}{2\pi}\ln\frac{r_s^2-d_i^2}{r_s r_{wi}}$, and $l_{ij}$ from the formula with $\cos(\theta_i-\theta_j)$ |

`GO TO (100,200,300) REFTYPE` is a computed GOTO that selects the case. The program then inverts **L** and sets $\mathbf{C}=\mu\varepsilon\,\mathbf{L}^{-1}$, valid because the medium is homogeneous.

### 4.2 `PCB.FOR` — lands on a board, no ground plane

**Method of moments** with pulse expansion. Each land of width w is split into `NCDIV` sub-strips of width δ = w/NCDIV, each with a constant charge density. Potential is enforced at the centre of each sub-strip, giving a system **A·ρ = V** of size N·NCDIV.

- `PHI(D)` — potential at distance D from the centre of a uniformly charged strip in free space, i.e. the integral of $-\frac{1}{2\pi\varepsilon_0}\ln|x|$ over the strip. D = 0 is the self term.
- `PHIB(D)` — same, on the surface of a dielectric board of thickness t. Uses an **image series** with coefficient $K=(\varepsilon_r-1)/(\varepsilon_r+1)$, summed until the change is below 1e-5.
- The "FILL A1 … FILL A MATRIX" loops exploit the **Toeplitz** structure: coefficients depend only on the distance between sub-strips, so only the first row of each block is computed and the rest is copied.
- After inverting A, the sum of the elements of block (i, j) times δ gives the **generalized capacitance matrix** $C^{g}_{ij}$ (charge vs. potential relative to infinity).
- Reduction to the line matrix with land `IREF` as reference:

$$C_{ij} = C^g_{ij} - \frac{\left(\sum_k C^g_{kj}\right)\left(\sum_k C^g_{ik}\right)}{\sum_k\sum_m C^g_{km}}$$

This is done twice: without dielectric to get **C₀** and then **L = C₀⁻¹/v₀²**, and with dielectric to get **C**.

**Input `PCB.IN`:** `N` (number of lands, including the reference), `NCDIV`, `IREF`, `W`, `SEP` (edge-to-edge separation), `T` (board thickness) — all in mils — and `ER`.
**Output:** L, C and C0.

### 4.3 `MSTRP.FOR` — coupled microstrip

Same structure as PCB (sub-strips, Toeplitz matrix, inversion), but with a **ground plane** under the substrate of thickness t:

- `PHI` includes the image of the sub-strip in the ground plane at depth 2t. Variables are normalized by δ/2: `DOW = 2D/δ`, `TOW = 4t/δ`.
- `PHIB` uses the image series of the dielectric slab over ground, with alternating terms $(-1)^n K^{n-1}$.
- The reference is the ground plane itself, so N lands give **N×N** matrices directly, with no reduction step. The sum of the blocks of $A^{-1}$ is already the C matrix.

**Input `MSTRP.IN`:** `N`, `NCDIV`, `W`, `SEP`, `T` (substrate thickness), `ER`.

### 4.4 `STRPLINE.FOR` — coupled stripline

Lands midway between two ground planes separated by `T`. The medium is **homogeneous**, so a single calculation suffices: **C = εr·C₀** and **L = C₀⁻¹/v₀²**.

`PHIB` numerically integrates the parallel-plate Green's function, $\ln\frac{u+1/u+2}{u+1/u-2}=\ln\coth^2\frac{\pi x}{2T}$ with $u=e^{\pi x/T}$, using **4-point Gauss–Legendre quadrature** (points ±0.339981 and ±0.861136, weights 0.652145 and 0.347855). The Gauss points never fall on x = 0, so the logarithmic singularity of the self term is avoided without special treatment.

**Input `STRPLINE.IN`:** `N`, `NCDIV`, `W`, `SEP`, `T` (separation between the planes), `ER`.

### 4.5 `RIBBON.FOR` — ribbon cable

The most numerically sophisticated program. Wires have radius `RW` and insulation thickness `TD`, with outer radius `RD = RW + TD`. The charge is expanded in a **Fourier series** around the perimeter with `NF` terms:

- free charge on the conductor surface (radius rw);
- bound polarization charge on the outer insulation surface (radius rd).

Because the wires lie in a single plane, the distribution is symmetric in θ and only **cosine** terms appear.

Boundary conditions are enforced at NF points on each perimeter. The angles are offset by a quarter step (`DELTA=ANGLE/4`) to avoid symmetry points.

- **A, B** — potential on the conductor surface due to free and bound charges.
- **C, D** — continuity of the normal component of **D** at the dielectric/air interface.

The coupled system is eliminated block-wise. Bound charges are $-\mathbf{D}^{-1}\mathbf{C}$ times the free charges, leading to $(\mathbf{A}-\mathbf{B}\mathbf{D}^{-1}\mathbf{C})^{-1}$. Only the zero-order coefficient carries net charge, so $C^g_{ij}=2\pi r_w\,a_0+2\pi r_d\,b_0$, and the code uses only row `(I-1)*NF+1`. The reduction with `IREF` is the same as in PCB.

**Input `RIBBON.IN`:** `N`, `NF`, `IREF`, `RW`, `TD`, `ER`, `S` (centre-to-centre separation). The program checks that the insulations do not overlap: 2(rw + td) ≤ s.

---

## 5. Stage 2 — SPICE models

Both read `SPICExxx.IN` with `NP1` (total number of conductors, N + 1) and `LEN` (line length in metres), and read **L** and **C** from `PUL.DAT`. External nodes are **101…10N** at the left end and **201…20N** at the right end; node 0 is the reference conductor.

### 5.1 `SPICEMTL.FOR` — exact modal model (lossless line)

This is the book's method for decoupling the MTL equations using real symmetric matrices:

1. `JACOBI` diagonalizes **C** with an orthogonal transformation U: $U^T C U = \theta^2$.
2. The symmetric matrix $\theta U^T L U \theta$ is formed and diagonalized with S.
3. The modal transformations **T_V** and **T_I** are obtained, with $T_I^{-1}=T_V^T$ and normalized, together with each mode's impedance $Z_{Cm}$ and velocity $v_m$. Subroutine `DIAG` returns these on the diagonals of C and L, which are destroyed.

The generated subcircuit contains:

- **V** — 0 V sources in series with the terminals, acting as ammeters.
- **EC…** with `POLY(N)` — controlled voltage sources, $\mathbf{V} = \mathbf{T}_V\,\mathbf{V}_m$, coefficients `TV(I,J)`.
- **FC…** with `POLY(N)` — controlled current sources, $\mathbf{I}_m=\mathbf{T}_V^{T}\mathbf{I}$, coefficients `TV(J,I)`.
- **T…** — N ideal, decoupled transmission lines, one per mode, with `Z0` = $Z_{Cm}$ and `TD` = LEN/$v_m$.

`JACOBI` is the threshold Jacobi method: it annihilates off-diagonal elements above the threshold, divides the threshold by N, and stops when it falls below 1e-6 of the initial norm. This is where the **arithmetic IFs** (`IF(x) 100,110,110`) appear; they branch on x being negative, zero or positive.

### 5.2 `SPICELPI.FOR` — lumped-Pi model

Valid only when the line is **electrically short** (length ≪ λ). For each conductor i:

- inductor `L` = $l_{ii}\cdot\mathcal{L}$ between nodes 10i and 20i;
- inductive coupling `K` = $l_{ij}/\sqrt{l_{ii}l_{jj}}$;
- at each end, a capacitor to ground equal to half of $\left(\sum_j c_{ij}\right)\mathcal{L}$. The off-diagonal $c_{ij}$ are negative, so this sum is the self capacitance to the reference;
- mutual capacitors `CM` = $-c_{ij}\mathcal{L}/2$ between conductors, at each end.

For a long line, cascade several Pi sections or use SPICEMTL.

---

## 6. Known issues and caveats

1. **`MSIZE=2`** in WIDESEP, SPICEMTL and SPICELPI: only 3 conductors are supported. For more, change the `PARAMETER` and recompile. The comment says `MSIZE=9`, but the code has 2. **SPICELPI does not check** `N > MSIZE` and simply overruns the array.
2. **RIBBON.FOR, FORMAT 44:** `'=C0(',I3,','I3,')'` is missing the comma before the second `I3`. `gfortran -std=legacy` should accept it with a warning; if it errors, change it to `',',I3,`.
3. **WIDESEP, case 3:** `R0=R0/CMTM` appears twice in the printing section. It does not affect results because R0 is not used afterwards.
4. In PCB, MSTRP and STRPLINE, a larger NCDIV improves accuracy, but N·NCDIV ≤ 100 (`MSIZE`), and inversion cost grows with the cube of that number.
5. The convergence tests of the series in `PHIB` are absolute (1e-5), not relative. With very high εr (K → 1) the series converges slowly.
6. The build instructions were derived from reading the source code; they have not been run through `gfortran` here. If a program fails to compile, check the compiler message against items 1–2 above first.

---

## 7. Publishing to GitHub

These steps assume the project lives in `~/emc` inside WSL, with the programs in `~/emc/cap9` and this file saved as `~/emc/cap9/README.md` (GitHub renders a file named `README.md` on the repository page).

### 7.1 Install and configure Git (once)

```bash
sudo apt install -y git
git config --global user.name "Your Name"
git config --global user.email "you@example.com"   # the e-mail linked to your GitHub account
git config --global init.defaultBranch main
```

### 7.2 Add a `.gitignore`

Keep compiled executables and per-run outputs out of the repository:

```bash
cd ~/emc
cat > .gitignore << 'EOF'
# Linux executables built from the .FOR sources
cap9/widesep
cap9/pcb
cap9/mstrp
cap9/ribbon
cap9/strpline
cap9/spicemtl
cap9/spicelpi
*.o
*.mod

# Outputs regenerated on every run
cap9/PUL.DAT
cap9/ref/

# Original DOS/Windows binaries
*.EXE

# VS Code (keep only tasks.json)
.vscode/*
!.vscode/tasks.json
EOF
```

Remove the corresponding lines if you want to keep the original `PUL.DAT` or the `.EXE` files in the repository.

### 7.3 Create the local repository and the first commit

```bash
cd ~/emc
git init
git add .
git status        # check what will be included
git commit -m "Crosstalk Fortran programs (Paul, ch. 9) with build setup and README"
```

### 7.4 Create the GitHub repository and push

**Option A — GitHub CLI (simplest):**

```bash
sudo apt install -y gh
gh auth login          # choose GitHub.com -> HTTPS -> log in with the browser
gh repo create emc --public --source=. --remote=origin --push
```

**Option B — website + SSH key:**

1. Generate a key and print its public part:
   ```bash
   ssh-keygen -t ed25519 -C "you@example.com"
   cat ~/.ssh/id_ed25519.pub
   ```
2. On GitHub, go to **Settings → SSH and GPG keys → New SSH key** and paste the key.
3. Create an empty **public** repository named `emc` at **github.com/new**, without a README, `.gitignore` or license, to avoid conflicts.
4. Connect and push:
   ```bash
   git remote add origin git@github.com:<your_user>/emc.git
   git push -u origin main
   ```

### 7.5 Future updates

From the WSL terminal:

```bash
git add .
git commit -m "Describe the change"
git push
```

In VS Code, the **Source Control** view (Ctrl+Shift+G) does the same graphically, using the Git installed in WSL.

### 7.6 Cloning on another machine

```bash
git clone https://github.com/<your_user>/emc.git
cd emc/cap9
make
```

### 7.7 Checking your GitHub authentication

**Logged-in accounts — GitHub CLI:**

```bash
gh auth status
```

Shows the authenticated accounts, which one is active, the Git protocol (HTTPS or SSH) and the token scopes. With more than one account, switch with `gh auth switch`.

**SSH key:**

```bash
ssh -T git@github.com
```

Replies with *"Hi \<user\>! You've successfully authenticated…"*, naming the account linked to the key.

**Commit identity:**

```bash
git config --global user.name
git config --global user.email
git config --list --show-origin | grep user   # also shows which file each value comes from
```

**Where `push` goes** (inside the project folder):

```bash
git remote -v
```

### 7.8 `gh` vs. SSH

Both are ways of **proving to GitHub who you are** when Git pushes and pulls. You only need one of them.

| | `gh` (HTTPS + token) | SSH (key pair) |
|---|---|---|
| How it authenticates | OAuth token (`gho_…`) stored in `~/.config/gh/hosts.yml` | Key pair `~/.ssh/id_ed25519`; the public key is registered on GitHub |
| Remote URL | `https://github.com/<your_user>/emc.git` | `git@github.com:<your_user>/emc.git` |
| What else it does | Also a GitHub command-line tool: create repositories, issues, PRs, releases, inspect Actions | Transports Git only; no access to the GitHub API |
| Password in daily use | None; the token is used automatically | The key passphrase, unless `ssh-agent` holds the key |

**Which one Git uses is decided by the remote URL**, not by preference: `https://…` uses the `gh` token, `git@github.com:…` uses the SSH key. A repository created with `gh repo create` gets the protocol configured in `gh` (HTTPS by default).

Switching between them:

```bash
git remote set-url origin git@github.com:<your_user>/emc.git     # use SSH
git remote set-url origin https://github.com/<your_user>/emc.git # back to HTTPS
```

To type the SSH passphrase only once per session:

```bash
eval "$(ssh-agent -s)" && ssh-add ~/.ssh/id_ed25519
```

**`user.name` and `user.email` are not authentication.** They are only the signature recorded in each commit:

1. GitHub links commits to your profile (avatar, contribution count) only if that e-mail is registered in **Settings → Emails**.
2. In a public repository the e-mail is visible in the history. To keep it private, use the *noreply* address shown on that same page:
   ```bash
   git config --global user.email "ID+<your_user>@users.noreply.github.com"
   ```

---

## 8. Copyright notice

The Fortran programs were written by Clayton R. Paul and are distributed with *Introduction to Electromagnetic Compatibility*, 3rd ed. (C. R. Paul, R. C. Scully, M. A. Steffka, Wiley, 2023). The `GAUSSJ` subroutine included in several programs is Copyright © 1987–1992 Numerical Recipes Software, whose notice grants permission for use only with the programs contained in the book. This repository adds build instructions and documentation; all rights to the original code remain with their respective holders.
