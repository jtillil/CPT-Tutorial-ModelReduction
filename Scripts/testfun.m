function out = testfun(in1, in2)

out = in1;

if exist("in2", "var")
    class(in2)
end

end