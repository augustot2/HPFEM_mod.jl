
function lagrange_aprox_spec(M,Q,nel,fun,a=-1,b=1,idir=[1,nel],nξ = 101)
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
    for e = 1:nel
        uu = fun(x[e])
        err = maxabs(uu-Ue[:,e])
        if err > maxerr maxerr = err end
    end

    return [maxerr HPFEM.nslvmodes(dof)]
    #return maxerr
end


function lagrange_oed_spec(M,Q,nel,fun,resp,a=-1,b=1,idir=[1 nel],nξ = 101)
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
  lixo = rand(1:72)
  cores  = gen_color()
  for e = 1:nel
      uu = resp(x[e])
      err = maxabs(uu-Ue[:,e])
      if err > maxerr maxerr = err end

  end

  return [maxerr HPFEM.nslvmodes(dof)]
  #return maxerr

end
