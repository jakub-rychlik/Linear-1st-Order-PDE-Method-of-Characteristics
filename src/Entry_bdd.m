% Generates points on entry part of boundary, where IVC is defined.
function entry_bdd = Entry_bdd(t)
   global Act_en_bdd;
   format long
   entry_bdd = [Act_en_bdd.xfunc(t);Act_en_bdd.yfunc(t)];
end 