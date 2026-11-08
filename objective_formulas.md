# Challenge objective formulas

## Objective

$$
J = \frac{E(\mathbf{p}_{\mathrm{edit}})}{E(\mathbf{p}_{\mathrm{rep}})}
+ C_{\mathrm{genes}}.
$$

## Recovery error

$$
E(\mathbf{p}) = \frac{E_T(\mathbf{p})^2 + E_\phi(\mathbf{p})^2}{2},
$$

$$
E_T(\mathbf{p}) = \frac{|T(\mathbf{p})-T_{\mathrm{WT}}|}{T_{\mathrm{WT}}},
\qquad
E_\phi(\mathbf{p}) = \sqrt{\frac{1}{5}\sum_{i=1}^{5}\delta_i(\mathbf{p})^2}.
$$

Here, each phase difference $\delta_i$ is measured in cycles, wrapped, and calculated after applying one common phase alignment to the five selected proteins.

## Gene-edit cost

$$
\begin{aligned}
C_{\mathrm{genes}} ={}& 0.405 I_{\mathrm{PER}}
+ 0.180 I_{\mathrm{CRY1}}
+ 0.156 I_{\mathrm{CRY2}} \\
&+ 0.168 I_{\mathrm{ROR}}
+ 0.192 I_{\mathrm{REV}}
+ 0.204 I_{\mathrm{BMAL}}.
\end{aligned}
$$

$$
I_g =
\begin{cases}
1, & \text{if any parameter in gene group }g\text{ differs from its reporter value},\\
0, & \text{otherwise}.
\end{cases}
$$

## Fully expanded objective

$$
\begin{aligned}
J ={}&
\frac{
\left(\dfrac{|T_{\mathrm{edit}}-T_{\mathrm{WT}}|}{T_{\mathrm{WT}}}\right)^2
+ \dfrac{1}{5}\sum_{i=1}^{5}\delta_{i,\mathrm{edit}}^2
}{
\left(\dfrac{|T_{\mathrm{rep}}-T_{\mathrm{WT}}|}{T_{\mathrm{WT}}}\right)^2
+ \dfrac{1}{5}\sum_{i=1}^{5}\delta_{i,\mathrm{rep}}^2
} \\
&+ 0.405 I_{\mathrm{PER}}
+ 0.180 I_{\mathrm{CRY1}}
+ 0.156 I_{\mathrm{CRY2}} \\
&+ 0.168 I_{\mathrm{ROR}}
+ 0.192 I_{\mathrm{REV}}
+ 0.204 I_{\mathrm{BMAL}}.
\end{aligned}
$$

The factors of $1/2$ in the recovery errors cancel in their ratio. The denominator is fixed by the reporter baseline.

## Multivariable function and gradient

Let $\mathbf{q}=(q_1,\ldots,q_k)$ contain the editable parameters being varied. Construct $\mathbf{p}(\mathbf{q})$ by inserting these values into the reporter vector and leaving all other entries fixed. Define

$$
f(\mathbf{q}) = J\bigl(\mathbf{p}(\mathbf{q})\bigr).
$$

For brevity, write $T(\mathbf{q})$ for the resulting period and $\delta_i(\mathbf{q})$ for the aligned, wrapped phase differences. Define the fixed denominator

$$
D =
\left(\frac{T_{\mathrm{rep}}-T_{\mathrm{WT}}}{T_{\mathrm{WT}}}\right)^2
+ \frac{1}{5}\sum_{i=1}^{5}\delta_{i,\mathrm{rep}}^2
= 2E(\mathbf{p}_{\mathrm{rep}}).
$$

Then

$$
f(\mathbf{q}) = \frac{1}{D}
\left[
\frac{(T(\mathbf{q})-T_{\mathrm{WT}})^2}{T_{\mathrm{WT}}^2}
+ \frac{1}{5}\sum_{i=1}^{5}\delta_i(\mathbf{q})^2
\right]
+ C_{\mathrm{genes}}(\mathbf{q}).
$$

Where the period and aligned phases are differentiable and the charged gene set is locally unchanged, the cost is constant. The chain rule gives

$$
\boxed{
\frac{\partial f}{\partial q_j}
= \frac{1}{D}
\left[
\frac{2(T(\mathbf{q})-T_{\mathrm{WT}})}{T_{\mathrm{WT}}^2}
\frac{\partial T}{\partial q_j}
+ \frac{2}{5}\sum_{i=1}^{5}
\delta_i(\mathbf{q})\frac{\partial\delta_i}{\partial q_j}
\right].
}
$$

Thus the gradient is the vector

$$
\nabla f(\mathbf{q}) =
\begin{bmatrix}
\dfrac{\partial f}{\partial q_1}\\
\vdots\\
\dfrac{\partial f}{\partial q_k}
\end{bmatrix}.
$$

Equivalently, under the same assumptions,

$$
\nabla f(\mathbf{q})
= \frac{1}{D}
\left[
\frac{2(T(\mathbf{q})-T_{\mathrm{WT}})}{T_{\mathrm{WT}}^2}\nabla T(\mathbf{q})
+ \frac{2}{5}\sum_{i=1}^{5}\delta_i(\mathbf{q})\nabla\delta_i(\mathbf{q})
\right].
$$

### Important limitations

- This is a symbolic chain-rule expression, not a numerical gradient. The period and phase sensitivities depend on the simulated model; the assignment asks you to implement finite-difference estimates.
- The gene-edit penalty is locally constant away from points where a gene becomes edited or unedited. At such switching points it jumps, so the full objective generally has no ordinary gradient. In particular, the reporter starting vector lies at such a switch for any gene you begin editing.
- Peak detection, phase wrapping, or loss of oscillations can also prevent smoothness. The expression above applies only where the relevant measurements are differentiable and the score is defined.

## Completion threshold

$$
J \leq 0.50,
\qquad
\text{Improvement (\%)} = 100(1-J).
$$
