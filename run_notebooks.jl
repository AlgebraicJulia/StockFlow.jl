function modinclude(filename)
    modname = gensym()
    @eval module $modname
      include($filename)
    end
end

function main()
  println("Starting to run notebooks ...")
  ipynb_files = filter(contains(r".jl$"), readdir("./jlexamples"; join=true))
  println("Checking: $ipynb_files")
  exit_code = 0
  errors = []
  ok = []
  for f in ipynb_files
      try
        modinclude(f)
      catch e
        push!(errors, (f, sprint(showerror, e)))
        exit_code = 1
      else
        push!(ok, f)
      end
  end
  println("OK: $ok")
  for (f, err_msg) in errors
      println("FAILED: $f")
      println("Reason: $err_msg")
  end
  exit(exit_code)
end

main()
