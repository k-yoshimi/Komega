Algorithm
=========

This library provides the four kinds of numerical solvers.
The kind of solvers is selected under the condition whether the Hamiltonian
:math:`{\hat H}` and/or the frequency :math:`z` are complex or real number.
It is noted that :math:`{\hat H}` must be Hermitian (symmetric)
for complex (real) number.

-  (:math:`{\hat H}`, :math:`z` ) = (complex, complex):
   Shifted Bi-Conjugate Gradient(BiCG) method :ref:`[1] <ref>`

-  (:math:`{\hat H}`, :math:`z` ) = (real, complex):
   Shifted Conjugate Orthogonal Conjugate Gradient(COCG) method :ref:`[2] <ref>`

-  (:math:`{\hat H}`, :math:`z` ) = (complex, real):
   Shifted Conjugate Gradient(CG) method (using complex vector)

-  (:math:`{\hat H}`, :math:`z` ) = (real, real):
   Shifted Conjugate Gradient(CG) method (using real vector)

For above methods, seed switching :ref:`[2] <ref>` is adopted.
Hereafter, the number of the left (right) side vector is
written as :math:`N_L` (:math:`N_R`).

The following methods can be selected by explicitly specifying them in "cg" area of the input file.

- QMR_SYM method :ref:`[3] <ref>`

- QMR_SYM(B) method :ref:`[3] <ref>`

The details of each algorithm are written as follows.

Shifted BiCG method with seed switching technique
-------------------------------------------------

:math:`G_{i j}(z_k) = 0 (i=1 \cdots N_L,\; j = 1 \cdots N_R,\; k=1 \cdots N_z)`

do :math:`j = 1 \cdots N_R`

   :math:`{\boldsymbol r} = {\boldsymbol \varphi_j}`,

   :math:`{\tilde {\boldsymbol r}} =` an arbitrary vector,
   :math:`{\boldsymbol r}^{\rm old} = {\tilde {\boldsymbol r}}^{\rm old} = {\bf 0}`

   :math:`p_{i k} = 0(i=1 \cdots N_L,\; k=1 \cdots N_z),\; \pi_k=\pi_k^{\rm old} = 1(k=1 \cdots N_z)`

   :math:`\rho = \infty,\; \alpha = 1,\; z_{\rm seed}=0`

   do iteration

      :math:`\circ` Seed equation

      :math:`\rho^{\rm old} = \rho,\; \rho = {\tilde {\boldsymbol r}}^* \cdot {\boldsymbol r}`

      :math:`\beta = \rho / \rho^{\rm old}`

      :math:`{\boldsymbol q} = (z_{\rm seed} {\hat I} - {\hat H}){\boldsymbol r}`

      :math:`\alpha^{\rm old} = \alpha,\; \alpha = \frac{\rho}{{\tilde {\boldsymbol r}}^*\cdot{\boldsymbol q} - \beta \rho / \alpha }`

      :math:`\circ` Shifted equation

      do :math:`k = 1 \cdots N_z`

         :math:`\pi_k^{\rm new} = [1+\alpha(z_k-z_{\rm seed})]\pi_k - \frac{\alpha \beta}{\alpha^{\rm old}}(\pi_k^{\rm old} - \pi_k)`

         do :math:`i = 1 \cdots N_L`

            :math:`p_{i k} = \frac{1}{\pi_k} {\boldsymbol \varphi}_i^* \cdot {\boldsymbol r} + \frac{\pi^{\rm old}_k \pi^{\rm old}_k}{\pi_k \pi_k} \beta p_{i k}`

            :math:`G_{i j}(z_k) = G_{i j}(z_k) + \frac{\pi_k}{\pi_k^{\rm new}} \alpha p_{i k}`

            :math:`\pi_k^{\rm old} = \pi_k`, :math:`\pi_k = \pi_k^{\rm new}`

         end do :math:`i`

      end do :math:`k`

      :math:`{\boldsymbol q} = \left( 1 + \frac{\alpha \beta}{\alpha^{\rm old}} \right) {\boldsymbol r} - \alpha {\boldsymbol q} - \frac{\alpha \beta}{\alpha^{\rm old}} {\boldsymbol r}^{\rm old},\; {\boldsymbol r}^{\rm old} = {\boldsymbol r},\; {\boldsymbol r} = {\boldsymbol q}`

      :math:`{\boldsymbol q} = (z_{\rm seed}^* {\hat I} - {\hat H}) {\tilde {\boldsymbol r}},\; {\boldsymbol q} = \left( 1 + \frac{\alpha^* \beta^*}{\alpha^{{\rm old}*}} \right) {\tilde {\boldsymbol r}} - \alpha^* {\boldsymbol q} - \frac{\alpha^* \beta^*}{\alpha^{{\rm old} *}} {\tilde {\boldsymbol r}}^{\rm old},\; {\tilde {\boldsymbol r}}^{\rm old} = {\tilde {\boldsymbol r}},\; {\tilde {\boldsymbol r}} = {\boldsymbol q}`

      :math:`\circ` Seed switch

      Search :math:`k` which gives the smallest :math:`|\pi_k|` . :math:`\rightarrow z_{\rm seed},\; \pi_{\rm seed},\; \pi_{\rm seed}^{\rm old}`

      :math:`{\boldsymbol r} = {\boldsymbol r} / \pi_{\rm seed},\; {\boldsymbol r}^{\rm old} = {\boldsymbol r}^{\rm old} / \pi_{\rm seed}^{\rm old},\; {\tilde {\boldsymbol r}} = {\tilde {\boldsymbol r}} / \pi_{\rm seed}^*,\; {\tilde {\boldsymbol r}}^{\rm old} = {\tilde {\boldsymbol r}}^{\rm old} / \pi_{\rm seed}^{{\rm old}*}`

      :math:`\alpha = (\pi_{\rm seed}^{\rm old} / \pi_{\rm seed}) \alpha`, :math:`\rho = \rho / (\pi_{\rm seed}^{\rm old} \pi_{\rm seed}^{\rm old})`

      :math:`\{\pi_k = \pi_k / \pi_{\rm seed},\; \pi_k^{\rm old} = \pi_k^{\rm old} / \pi_{\rm seed}^{\rm old}\}`

      if( :math:`|{\boldsymbol r}| <` Threshold) exit

   end do iteration

end do :math:`j`

Shifted COCG method with seed switching technique
-------------------------------------------------

This method is obtained by
:math:`{\tilde {\boldsymbol r}} = {\boldsymbol r}^*,\; {\tilde {\boldsymbol r}}^{\rm old} = {\boldsymbol r}^{{\rm old}*}`
in the BiCG method.

:math:`G_{i j}(z_k) = 0 (i=1 \cdots N_L,\; j = 1 \cdots N_R,\; k=1 \cdots N_z)`

do :math:`j = 1 \cdots N_R`

   :math:`{\boldsymbol r} = {\boldsymbol \varphi_j}`, :math:`{\boldsymbol r}^{\rm old} = {\bf 0}`

   :math:`p_{i k} = 0(i=1 \cdots N_L,\; k=1 \cdots N_z),\; \pi_k=\pi_k^{\rm old} = 1(k=1 \cdots N_z)`

   :math:`\rho = \infty,\; \alpha = 1,\; z_{\rm seed}=0`

   do iteration

      :math:`\circ` Seed equation

      :math:`\rho^{\rm old} = \rho,\; \rho = {\boldsymbol r} \cdot {\boldsymbol r}`

      :math:`\beta = \rho / \rho^{\rm old}`

      :math:`{\boldsymbol q} = (z_{\rm seed} {\hat I} - {\hat H}){\boldsymbol r}`

      :math:`\alpha^{\rm old} = \alpha,\; \alpha = \frac{\rho}{{\boldsymbol r}\cdot{\boldsymbol q} - \beta \rho / \alpha }`

      :math:`\circ` Shifted equation

      do :math:`k = 1 \cdots N_z`

         :math:`\pi_k^{\rm new} = [1+\alpha(z_k-z_{\rm seed})]\pi_k - \frac{\alpha \beta}{\alpha^{\rm old}}(\pi_k^{\rm old} - \pi_k)`

         do :math:`i = 1 \cdots N_L`

            :math:`p_{i k} = \frac{1}{\pi_k} {\boldsymbol \varphi}_i^* \cdot {\boldsymbol r} + \frac{\pi^{\rm old}_k \pi^{\rm old}_k}{\pi_k \pi_k} \beta p_{i k}`

            :math:`G_{i j}(z_k) = G_{i j}(z_k) + \frac{\pi_k}{\pi_k^{\rm new}} \alpha p_{i k}`

            :math:`\pi_k^{\rm old} = \pi_k`, :math:`\pi_k = \pi_k^{\rm new}`

         end do :math:`i`

      end do :math:`k`

      :math:`{\boldsymbol q} = \left( 1 + \frac{\alpha \beta}{\alpha^{\rm old}} \right) {\boldsymbol r} - \alpha {\boldsymbol q} - \frac{\alpha \beta}{\alpha^{\rm old}} {\boldsymbol r}^{\rm old},\; {\boldsymbol r}^{\rm old} = {\boldsymbol r},\; {\boldsymbol r} = {\boldsymbol q}`

      :math:`\circ` Seed switch

      Search :math:`k` which gives the smallest :math:`|\pi_k|` . :math:`\rightarrow z_{\rm seed},\; \pi_{\rm seed},\; \pi_{\rm seed}^{\rm old}`
                  
      :math:`{\boldsymbol r} = {\boldsymbol r} / \pi_{\rm seed},\; {\boldsymbol r}^{\rm old} = {\boldsymbol r}^{\rm old} / \pi_{\rm seed}^{\rm old}`

      :math:`\alpha = (\pi_{\rm seed}^{\rm old} / \pi_{\rm seed}) \alpha`, :math:`\rho = \rho / (\pi_{\rm seed}^{\rm old} \pi_{\rm seed}^{\rm old})`

      :math:`\{\pi_k = \pi_k/\pi_{\rm seed},\; \pi_k^{\rm old} = \pi_k^{\rm old} / \pi_{\rm seed}^{\rm old}\}`

      if( :math:`|{\boldsymbol r}| <` Threshold) exit

   end do iteration

end do :math:`j`

Shifted CG method with seed switching technique
-----------------------------------------------

This method is obtained by
:math:`{\tilde {\boldsymbol r}} = {\boldsymbol r},\; {\tilde {\boldsymbol r}}^{\rm old} = {\boldsymbol r}^{\rm old}`
in the BiCG method.

:math:`G_{i j}(z_k) = 0 (i=1 \cdots N_L,\; j = 1 \cdots N_R,\; k=1 \cdots N_z)`

do :math:`j = 1 \cdots N_R`

   :math:`{\boldsymbol r} = {\boldsymbol \varphi_j}`, :math:`{\boldsymbol r}^{\rm old} = {\bf 0}`

   :math:`p_{i k} = 0(i=1 \cdots N_L,\; k=1 \cdots N_z),\; \pi_k=\pi_k^{\rm old} = 1(k=1 \cdots N_z)`

   :math:`\rho = \infty,\; \alpha = 1,\; z_{\rm seed}=0`

   do iteration

      :math:`\circ` Seed equation

      :math:`\rho^{\rm old} = \rho,\; \rho = {\boldsymbol r}^* \cdot {\boldsymbol r}`

      :math:`\beta = \rho / \rho^{\rm old}`

      :math:`{\boldsymbol q} = (z_{\rm seed} {\hat I} - {\hat H}){\boldsymbol r}`

      :math:`\alpha^{\rm old} = \alpha,\; \alpha = \frac{\rho}{{\boldsymbol r}^* \cdot {\boldsymbol q} - \beta \rho / \alpha }`

      :math:`\circ` Shifted equation

      do :math:`k = 1 \cdots N_z`

         :math:`\pi_k^{\rm new} = [1+\alpha(z_k-z_{\rm seed})]\pi_k - \frac{\alpha \beta}{\alpha^{\rm old}}(\pi_k^{\rm old} - \pi_k)`

         do :math:`i = 1 \cdots N_L`

            :math:`p_{i k} = \frac{1}{\pi_k} {\boldsymbol \varphi}_i^* \cdot {\boldsymbol r} + \left(\frac{\pi^{\rm old}_k}{\pi_k } \right)^2 \beta p_{i k}`

            :math:`G_{i j}(z_k) = G_{i j}(z_k) + \frac{\pi_k}{\pi_k^{\rm new}} \alpha p_{i k}`

            :math:`\pi_k^{\rm old} = \pi_k`, :math:`\pi_k = \pi_k^{\rm new}`

         end do :math:`i`

      end do :math:`k`

      :math:`{\boldsymbol q} = \left( 1 + \frac{\alpha \beta}{\alpha^{\rm old}} \right) {\boldsymbol r} - \alpha {\boldsymbol q} - \frac{\alpha \beta}{\alpha^{\rm old}} {\boldsymbol r}^{\rm old},\; {\boldsymbol r}^{\rm old} = {\boldsymbol r},\; {\boldsymbol r} = {\boldsymbol q}`

      :math:`\circ` Seed switch

      Search :math:`k` which gives the minimum value of :math:`|\pi_k|` . :math:`\rightarrow z_{\rm seed},\; \pi_{\rm seed},\; \pi_{\rm seed}^{\rm old}`

      :math:`{\boldsymbol r} = {\boldsymbol r} / \pi_{\rm seed},\; {\boldsymbol r}^{\rm old} = {\boldsymbol r}^{\rm old} / \pi_{\rm seed}^{\rm old}`

      :math:`\alpha = (\pi_{\rm seed}^{\rm old} / \pi_{\rm seed}) \alpha`, :math:`\rho = \rho / {\pi_{\rm seed}^{\rm old}}^2`

      :math:`\{\pi_k = \pi_k/\pi_{\rm seed},\; \pi_k^{\rm old} = \pi_k^{\rm old}/\pi_{\rm seed}^{\rm old}\}`

      if( :math:`|{\boldsymbol r}| <` Threshold) exit

   end do iteration

end do :math:`j`

QMR_SYM method
------------------

:math:`\boldsymbol{x}_{0}^{(\ell)}=\boldsymbol{p}_{-1}^{(\ell)}=\boldsymbol{p}_{0}^{(\ell)}=0, \boldsymbol{v}_{1}=\boldsymbol{b}/(\boldsymbol{b}^{T}\boldsymbol{b})^{1/2}, g_{1}^{(\ell)}=(\boldsymbol{b}^{T}\boldsymbol{b})^{1/2}` 

do :math:`n = 1, 2, \cdots`

   :math:`\circ` The complex symmetric Lanczos process

   :math:`\alpha_n = \boldsymbol{v}_{n}^{T}A\boldsymbol{v}_{n}`

   :math:`\tilde{\boldsymbol{v}}_{n+1}=A\boldsymbol{v}_n-\alpha_{n}\boldsymbol{v}_n-\beta_{n-1}\boldsymbol{v}_{n-1}`

   :math:`\beta_{n}=(\tilde{\boldsymbol{v}}^T_{n+1}\tilde{\boldsymbol{v}}_{n+1})^{1/2}`

   :math:`\boldsymbol{v}_{n+1}=\tilde{\boldsymbol{v}}_{n+1}/\beta_{n}`

   :math:`t^{(\ell)}_{n-1,n}=\beta_{n-1}, t^{(\ell)}_{n,n}=\alpha_{n}+\sigma_{\ell}, t^{(\ell)}_{n+1,n}=\beta_{n}`

   :math:`\circ` Solve least squares problems by Givens rotations

   do :math:`\ell = 1, 2, \cdots, m`

      if( :math:`||\boldsymbol{r}^{(\ell)}_n||_2/||\boldsymbol{b}||_2\geq\epsilon` )

         do :math:`i=\rm{max}\{1,n-2\},\cdot,n-1`

            :math:`\left[\begin{array}{c}{t_{i, n}^{(\ell)}} \\{t_{i+1, n}^{(\ell)}}\end{array}\right]=\left[\begin{array}{cc}{c_{i}^{(\ell)}} & {s_{i}^{(\ell)}} \\{-\bar{s}_{i}^{(\ell)}} & {c_{i}^{(\ell)}}\end{array}\right]\left[\begin{array}{c}{t_{i, n}^{(\ell)}} \\{t_{i+1}^{(\ell)}}\end{array}\right]`

         end do :math:`i`

         :math:`c_{n}^{(\ell)}=\frac{\left|t_{n, n}^{(\ell)}\right|}{\sqrt{\left|t_{n, n}^{(\ell)}\right|^{2}+\left|t_{n+1, n}^{(\ell)}\right|^{2}}}`

         :math:`\bar{s}_{n}^{(\ell)}=\frac{t_{n+1, n}^{(\ell)}}{t_{n, n}^{(\ell)}} c_{n}^{(\ell)}`

         :math:`t_{n, n}^{(\ell)}=c_{n}^{(\ell)} t_{n, n}^{(\ell)}+s_{n}^{(\ell)} t_{n+1, n}^{(\ell)}`

         :math:`\left[\begin{array}{c}{g_{n}^{(\ell)}} \\{g_{n+1}^{(\ell)}}\end{array}\right]=\left[\begin{array}{cc}{c_{n}^{(\ell)}} & {s_{n}^{(\ell)}} \\{-\bar{s}_{n}^{(\ell)}} & {c_{n}^{(\ell)}}\end{array}\right]\left[\begin{array}{c}{g_{n}^{(\ell)}} \\{0}\end{array}\right]`

         :math:`\circ` Update approximate solutions :math:`x_{n}^{(\ell)}`

         :math:`\boldsymbol{p}_{n}^{(\ell)}=\boldsymbol{v}_{n}-\left(t_{n-2, n}^{(\ell)} / t_{n-2, n-2}^{(\ell)}\right) \boldsymbol{p}_{n-2}^{(\ell)}-\left(t_{n-1, n}^{(\ell)} / t_{n-1, n-1}^{(\ell)}\right) \boldsymbol{p}_{n-1}^{(\ell)}`

         :math:`\boldsymbol{x}_{n}^{(\ell)}=\boldsymbol{x}_{n-1}^{(\ell)}+\left(g_{n}^{(\ell)} / t_{n, n}^{(\ell)}\right) \boldsymbol{p}_{n}^{(\ell)}`

      endif

   end do :math:`\ell`

   if( :math:`||\boldsymbol{r}^{(\ell)}_n||_2/||\boldsymbol{b}||_2\leq\epsilon` for all :math:`\ell` ) then exit.

end do :math:`n`

QMR_SYM(B) method
--------------------

:math:`\boldsymbol{x}_{0}^{(\ell)}=\boldsymbol{p}_{-1}^{(\ell)}=\boldsymbol{p}_{0}^{(\ell)}=0, \boldsymbol{v}_{1}=\boldsymbol{b}/(\boldsymbol{b}^{T}\boldsymbol{b})^{1/2}, g_{1}^{(\ell)}=(\boldsymbol{b}^{T}\boldsymbol{b})^{1/2}` 

do :math:`n = 1, 2, \cdots`

   :math:`\circ` The complex symmetric Lanczos process

   :math:`\alpha_n = \boldsymbol{v}_{n}^{T}A\boldsymbol{v}_{n}`

   :math:`\tilde{\boldsymbol{v}}_{n+1}=A\boldsymbol{v}_n-\alpha_{n}\boldsymbol{v}_n-\beta_{n-1}\boldsymbol{v}_{n-1}`

   :math:`\beta_{n}=(\tilde{\boldsymbol{v}}^T_{n+1}\tilde{\boldsymbol{v}}_{n+1})^{1/2}`

   :math:`\boldsymbol{v}_{n+1}=\tilde{\boldsymbol{v}}_{n+1}/\beta_{n}`

   :math:`t^{(\ell)}_{n-1,n}=\beta_{n-1}, t^{(\ell)}_{n,n}=\alpha_{n}+\sigma_{\ell}, t^{(\ell)}_{n+1,n}=\beta_{n}`

   :math:`\circ` Solve weighted least squares problems

   do :math:`\ell = 1, 2, \cdots, m`

      if( :math:`||\boldsymbol{r}^{(\ell)}_n||_2/||\boldsymbol{b}||_2\geq\epsilon` )

         do :math:`i=\rm{max}\{1,n-1\},\cdot,n-1`

            :math:`t_{i+1, n}^{(\ell)}=f_{i}^{(\ell)} t_{i, n}^{(\ell)}+t_{i+1, n}^{(\ell)}`

         end do :math:`i`

         :math:`f_{n}^{(\ell)}=-\frac{t_{n+1, n}^{(\ell)}}{t_{n, n}^{(\ell)}}`

         :math:`t_{n+1, n}^{(\ell)}=0`

         :math:`\widetilde{g}_{n+1}^{(\ell)}=f_{n}^{(\ell)} \widetilde{g}_{n}^{(\ell)}`

         :math:`\circ` Update approximate solutions :math:`x_{n}^{(\ell)}`

         :math:`p_{n}^{(\ell)}=v_{n}-\left(t_{n-1, n}^{(\ell)} / t_{n-1, n-1}^{(\ell)}\right) p_{n-1}^{(\ell)}`

         :math:`\boldsymbol{x}_{n}^{(\ell)}=\boldsymbol{x}_{n-1}^{(\ell)}+\left(\tilde{g}_{n}^{(\ell)} / t_{n, n}^{(\ell)}\right) \boldsymbol{p}_{n}^{(\ell)}`

      endif

   end do :math:`\ell`

   if( :math:`||\boldsymbol{r}^{(\ell)}_n||_2/||\boldsymbol{b}||_2\leq\epsilon` for all :math:`\ell` ) then exit.

end do :math:`n`

