    else
          thisname = class(val);  
       end
    end
    
    % Remove '.' characters from variable name 
    thisname = strrep(thisname,'.','_');
    
    % See if this variable type is already present in the 
    % variable list
    namelist = get(hVariableTable,'VariableNameList');
    namelistcount = get(hVariableTable,'VariableNameListCount'); 
    ind = find(strcmpi(namelist,thisname)==true);
    
    % If it is not in the list, then add it
    if isempty(ind) 
        count = 1;
        newname = sprintf('%s%d',thisname,count);
        set(hVariableTable,'VariableNameList',{namelist{:},thisname});
        set(hVariableTable,'VariableNameListCount',[namelistcount,count]);
    
    % If it is in the list, increment variable count
    else
        count = namelistcount(ind(1))+1;
        namelistcount(ind(1)) = count;
        set(hVariableTable,'VariableNameListCount',namelistcount);
        newname = sprintf('%s%d',thisname,count);
    end
    set(hArg,'String',newname);
end                    5d  h   6d  p   7d  x   8d     9d     :d     ;d     <d      =d  ¨   >d  °   ?d  ¸   @d  À   Ad  È   Bd  Ð   Cd  Ø   Dd  à   Ed  è   Fd  ð   Gd  ø   Hd      Id     Jd     Kd     Ld      Md  (   Nd  0   Od  8   Pd  @   Qd  H   Rd  P   Sd  X   Td  `   Ud  h   Vd  p   Wd  x   Xd     Yd     Zd     [d     \d      ]d  ¨   ^d  °   _d  ¸   `d  À   ad  È   bd  Ð   cd  Ø   dd  à   ed  è   fd  ð   gd  ø   hd      id     jd     kd     ld      md  (   nd  0   od  8   pd  @   qd  H   rd  P   sd  X   td  `   ud  h   vd  p   wd  x   xd     yd     zd     {d     |d      }d  ¨   ~d  °   d  ¸   d  À   d  È   d  Ð   d  Ø   d  à   d  è   d  ð   d  ø   d      d     d     d     d      d  (   d  0   d  8   d  @   d  H   d  P   d  X   d  `   d  h   d  p   d  x   d     d     d     d     d      d  ¨   d  °   d  ¸    d  À   ¡d  È   ¢d  Ð   £d  Ø   ¤d  à   ¥d  è   ¦d  ð   §d  ø   ¨d      ©d     ªd     «d     ¬d      ­d  (   ®d  0   ¯d  8   °d  @   ±d  H   ²d  P   ³d  X   ´d  `   µd  h   ¶d  p   ·d  x   ¸d     ¹d     ºd     »d     ¼d      ½d  ¨   ¾d  °   ¿d  ¸   Àd  À   Ád  È   Âd  Ð   Ãd  Ø   Äd  à   Åd  è   Æd  ð   Çd  ø   Èd      Éd     Êd     Ëd     Ìd      Íd  (   Îd  0   Ïd  8   Ðd  @   Ñd  H   Òd  P  