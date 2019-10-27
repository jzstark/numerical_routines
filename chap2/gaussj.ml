(* Linear equation solution by Gauss-Jordan elimination *)
let gaussj a b =
  (*TODO: shape check *)
  let a = M.copy a in
  let b = M.copy b in

  let n = (M.shape a).(0) in
  let m = (M.shape b).(1) in
  let icol = ref 0 in
  let irow = ref 0 in
  let dum = ref 0.0 in
  let pivinv = ref 0.0 in
  let indxc = Array.make n 0 in
  let indxr = Array.make n 0 in
  let ipiv  = Array.make n 0 in
  for i = 0 to n - 1  do
    let big = ref 0.0 in
    for j = 0 to n - 1 do
      if ipiv.(j) != 1 then (
        for k = 0 to n - 1 do
          if ipiv.(k) == 0 then (
            let v = M.get a [|j; k|] |> abs_float in
            if (v >= !big) then (
              big := v;
              irow := j; icol := k;
            )
          )
        done
      )
    done;
    ipiv.(!icol) <- ipiv.(!icol) + 1;

    (* TODO: how to swap elegantly? *)
    if (!irow != !icol) then (
      for l = 0 to n - 1 do
        let u = M.get a [|!irow; l|] in
        let v = M.get a [|!icol; l|] in
        M.set a [|!icol; l|] u;
        M.set a [|!irow; l|] v
      done;

      for l = 0 to n - 1 do
        let u = M.get b [|!irow; l|] in
        let v = M.get b [|!icol; l|] in
        M.set b [|!icol; l|] u;
        M.set b [|!irow; l|] v
      done
    );

    indxr.(i) <- !irow;
    indxc.(i) <- !icol;
    let p = M.get a [|!icol; !icol|] in
    if (p == 0.0) then failwith "gaussj: Singular Matrix";
    pivinv :=  1.0 /. p;
    M.set a [|!icol; !icol|] 1.0;
    for l = 0 to n - 1 do
      let prev = M.get a [|!icol; l|] in
      M.set a [|!icol; l|] (prev *. !pivinv)
    done;
    for l = 0 to m - 1 do
      let prev = M.get b [|!icol; l|] in
      M.set b [|!icol; l|] (prev *. !pivinv)
    done;

    for ll = 0 to n - 1 do
      if (ll != !icol) then (
        dum := M.get a [|ll; !icol|];
        M.set a [|ll; !icol|] 0.0;
        for l = 0 to n - 1 do
          let p = M.get a [|!icol; l|] in
          let prev = M.get a [|ll; l|] in
          M.set a [|ll; l|] (prev -. p *. !dum)
        done;
        for l = 0 to n - 1 do
          let p = M.get b [|!icol; l|] in
          let prev = M.get b [|ll; l|] in
          M.set b [|ll; l|] (prev -. p *. !dum)
        done
      )
    done

  done;

  for l = n - 1 downto 0 do
    if (indxr.(l) != indxc.(l)) then (
      for k = 0 to n - 1 do
        let u = M.get a [|k; indxr.(l)|] in
        let v = M.get a [|k; indxc.(l)|] in
        M.set a [|k; indxc.(l)|] u;
        M.set a [|k; indxr.(l)|] v
      done
    )
  done;

  a, b
