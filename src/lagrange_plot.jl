function gen_color()
    letter = num2hex(rand(0:255))[15:16]
    for i in 1:2
        letter = string(letter,num2hex(rand(0:255))[15:16])
    end
    cor = string("#",letter)
    return cor
end

function lagrange_oed_plot(M,Q,nel,fun,resp,a=-1,b=1,idir=[1 nel],nξ = 101)
  nnodes = nel + 1
  nodes = collect(linspace(a, b, nnodes));

  lagr = HPFEM.Lagrange1d(M);
  quad = HPFEM.QuadType(Q);
  base = HPFEM.Basis1d(lagr, quad);

  lmap = HPFEM.locmap(base)
  dof = HPFEM.DofMap1d(lmap, nnodes, idir);

  elemento = [HPFEM.Element1d(e, nodes[e], nodes[e+1], base) for e = 1:nel]
  solver = HPFEM.CholeskySC(dof, HPFEM.BBSymTri);

  for e = 1:nel
      x  = elemento[e].x

      Ae = HPFEM.mass_matrix(base, elemento[e])
      Se = HPFEM.stiff_matrix(base,elemento[e])
      Ae = Ae + Se
      HPFEM.add_local_matrix(solver, e, Ae)
  end

  bnd = HPFEM.bndidx(lmap)

  Fe = zeros(HPFEM.nmodes(lmap), nel)

  for e = 1:nel
      fe = fun(elemento[e].x)
      HPFEM.add_rhs!(base, elemento[e], fe, sub(Fe, :, e))
  end
  Fe[bnd[1],1] = resp(a)
  Fe[bnd[2],nel]= resp(b)

  HPFEM.solve!(solver, Fe)

  nξ = 101
  ξ = collect(linspace(-1,1,nξ));
  ϕ = zeros(nξ, M)
  for i = 1:M
      ϕ[:,i] = base(ξ, i)
  end

  Ue = ϕ * Fe

  x = [(1-ξ)*el.a/2 + (1+ξ)*el.b/2 for el in elemento]
  maxerr = -1000000

  cores  = gen_color()
  for e = 1:nel
      uu = resp(x[e])
      err = maxabs(uu-Ue[:,e])
      if err > maxerr maxerr = err end
      plot( x[e],uu, "-b")
      plot(x[e],Ue[:,e], "--",color = cores)
      #legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
  end
  suptitle("M = $M , Q=$Q, Nel = $nel")
  title("erro = $maxerr")
  return maxerr

end

function lagrange_aprox_plot(M,Q,nel,fun,a=-1,b=1,idir=[1,nel],nξ = 101)
    nnodes = nel + 1
    lagr = HPFEM.Lagrange1d(M);
    quad = HPFEM.QuadType(Q);
    base = HPFEM.Basis1d(lagr, quad);
    nodes = collect(linspace(a, b, nnodes));
    lmap = HPFEM.locmap(base)
    dof = HPFEM.DofMap1d(lmap, nnodes, idir);

    elemento = [HPFEM.Element1d(e, nodes[e], nodes[e+1], base) for e = 1:nel]
    solver = HPFEM.CholeskySC(dof, HPFEM.BBSymTri);

    for e = 1:nel
        Ae = HPFEM.mass_matrix(base, elemento[e])
        HPFEM.add_local_matrix(solver, e, Ae)
    end

    Fe = zeros(HPFEM.nmodes(lmap), nel)
    for e = 1:nel
        fe = fun(elemento[e].x)
        HPFEM.add_rhs!(base, elemento[e], fe, sub(Fe, :, e))
    end

    HPFEM.solve!(solver, Fe)

    nξ = 101
    ξ = collect(linspace(-1,1,nξ));
    ϕ = zeros(nξ, M)
    for i = 1:M
        ϕ[:,i] = base(ξ, i)
    end
    Ue = ϕ * Fe;
    x = [(1-ξ)*el.a/2 + (1+ξ)*el.b/2 for el in elemento]
    maxerr = -1000


    cores  = gen_color()
    for e = 1:nel
        uu = fun(x[e])
        err = maxabs(uu-Ue[:,e])
        if err > maxerr maxerr = err end
        plot( x[e],uu, "-b")
        plot(x[e],Ue[:,e], "--",color = cores)
        #legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    end
    #suptitle("M = $M , Q=$Q, Nel = $nel")
    title("N = $M , Nel = $nel")
    #title("erro = $maxerr")

    return [maxerr HPFEM.nslvmodes(dof)]

    return 0

end



function lagrange_oed_neum_plot(M,Q,nel,rhsfun,uexact,duexact,a=-1,b=1,idir=[1,nel],nξ = 101)
  TF = Float64
  nnodes = nel + 1
  b1 = HPFEM.Lagrange1d(M, TF)
  quad = HPFEM.QuadType(Q, HPFEM.GLJ, TF)
  bas = HPFEM.Basis1d(b1,quad, TF)
  lmap = HPFEM.locmap(bas)
  dof = HPFEM.DofMap1d(lmap, nnodes, idir);

  nodes = [TF(x) for x in linspace(a, b, nnodes)];

  elems = [HPFEM.Element1d(e, nodes[e], nodes[e+1], bas) for e = 1:nel];

  solver = HPFEM.CholeskySC(dof, HPFEM.BBMatrix1d, TF);



  for e = 1:nel
      Ae = zeros(TF, M, M)
      HPFEM.add_stiff_matrix!(bas, elems[e], Ae)
      HPFEM.add_mass_matrix!(bas, elems[e], Ae)
      HPFEM.add_local_matrix(solver, e, Ae)
  end


  solver.Abb

  Fe = zeros(HPFEM.nmodes(lmap), nel)
  bnd = HPFEM.bndidx(lmap)

  for e = 1:nel
      fe = rhsfun(elems[e].x)
      HPFEM.add_rhs!(bas, elems[e], fe, sub(Fe, :, e))
  end
  Fe[bnd[2],nel] += bas(one(TF), bnd[2]) * duexact(b) ### A: Adicionar o termo vdu/dx no lado direito!
  # Apply Dirichilet BCs:
  Fe[bnd[1],1] = uexact(a);

  HPFEM.solve!(solver, Fe);

  nξ = 101
  ξ = collect(linspace(-1,1,nξ));
  ϕ = zeros(nξ, M)
  for i = 1:M
      ϕ[:,i] = bas(ξ, i)
  end

  Ue = ϕ * Fe;

  #plot
  maxerr = 0.0
  for e = 1:nel
      el = elems[e]
      x = (1-ξ)*el.a/2 + (1+ξ)*el.b/2
      uu = uexact(x)
      err = maxabs(uu-Ue[:,e])
      if err > maxerr maxerr = err end

      plot(x, Ue[:,e], "r", x, uu, "b")
  end
  return maxerr

end
