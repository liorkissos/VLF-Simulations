 function Store_Data1(src,event)
 
 %disp('Inside Store_Data1')
    global BBB
    global CCC

     BBB = [BBB;event.Data]; % operational useful variable. rests at every session call
     CCC = [CCC;event.Data]; % testability variable. accumulates all the recorded sessions so far
    
     
  %   stop(src)

 end

 