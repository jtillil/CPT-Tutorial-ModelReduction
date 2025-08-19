function odejacstr = odejacsym2odejacstr(odejacsym, model)

I = model.I;
odejacstr = "";

for i = 1:I.nstates
for j = 1:I.nstates
% if j == 1
%     odejacstr = odejacstr + "%{ Variable " + string(i) + " %}\n";
% end
if ~isequal(odejacsym(i, j), sym(0))
    odejacstr = odejacstr + "DF(I." + I.nmstate{i} + ",I." + I.nmstate{j} + ") = ";
    currstr = string(odejacsym(i, j));
    for name = I.nmpar
        occurances = strfind(currstr,name{1});
        for pos = flip(occurances)
            if ~(pos+strlength(name{1}) > strlength(currstr)) && pos ~= 1
                if contains(extract(currstr, pos+strlength(name{1})), ["+", "-", "*", "/", "^", " ", "(", ")"]) && contains(extract(currstr, pos-1), ["+", "-", "*", "/", "^", " ", "(", ")"])
                    currstr = replaceBetween(currstr, pos, pos+strlength(name{1})-1, "par(I." + name{1} + ")");
                end
            elseif pos+strlength(name{1}) == strlength(currstr)+1 && pos ~= 1
                if contains(extract(currstr, pos-1), ["+", "-", "*", "/", "^", " ", "(", ")"])
                    currstr = replaceBetween(currstr, pos, pos+strlength(name{1})-1, "par(I." + name{1} + ")");
                end
            elseif pos+strlength(name{1}) == strlength(currstr)+1 && pos == 1
                currstr = replaceBetween(currstr, pos, pos+strlength(name{1})-1, "par(I." + name{1} + ")");
            else
                if contains(extract(currstr, pos+strlength(name{1})), ["+", "-", "*", "/", "^", " ", "(", ")"])
                    currstr = replaceBetween(currstr, pos, pos+strlength(name{1})-1, "par(I." + name{1} + ")");
                end
            end
        end
    end
    for name = I.nmstate
        occurances = strfind(currstr,name{1});
        for pos = flip(occurances)
            if ~(pos+strlength(name{1}) > strlength(currstr)) && pos ~= 1
                if contains(extract(currstr, pos+strlength(name{1})), ["+", "-", "*", "/", "^", " ", "(", ")"]) && contains(extract(currstr, pos-1), ["+", "-", "*", "/", "^", " ", "(", ")"])
                    currstr = replaceBetween(currstr, pos, pos+strlength(name{1})-1, "X(I." + name{1} + ")");
                end
            elseif pos+strlength(name{1}) == strlength(currstr)+1 && pos ~= 1
                if contains(extract(currstr, pos-1), ["+", "-", "*", "/", "^", " ", "(", ")"])
                    currstr = replaceBetween(currstr, pos, pos+strlength(name{1})-1, "X(I." + name{1} + ")");
                end
            elseif pos+strlength(name{1}) == strlength(currstr)+1 && pos == 1
                currstr = replaceBetween(currstr, pos, pos+strlength(name{1})-1, "X(I." + name{1} + ")");
            else
                if contains(extract(currstr, pos+strlength(name{1})), ["+", "-", "*", "/", "^", " ", "(", ")"])
                    currstr = replaceBetween(currstr, pos, pos+strlength(name{1})-1, "X(I." + name{1} + ")");
                end
            end
        end
    end
    odejacstr = odejacstr + currstr;
    odejacstr = odejacstr + ";\n";
end
if j == I.nstates
    odejacstr = odejacstr + "\n";
end
end
end

fprintf(odejacstr);

end