function modinclude(filename)
    modname = gensym()
    @eval module $modname
      include($filename)
    end
end

function main()
  println("Starting to run notebooks ...")
  exit_code = 0
  errors = []
  ok = []
  for (root, dirs, files) in walkdir("./jlexamples")
      if !isempty(files)
          println("Checking examples in $root")
          for f in files
              try
                  modinclude(root * "/" * f)
              catch e
                  push!(errors, (root * "/" * f, sprint(showerror, e)))
                  exit_code = 1
              else
                  push!(ok, f)
              end
          end
      end
  end
  println("OK: $ok")
  failed_files = [f for (f, _) in errors]
  println("FAILED: $failed_files")
  for (f, err_msg) in errors
      println("FAILED: $f")
      println("Reason: $err_msg")
  end
  exit(exit_code)
end

main()
