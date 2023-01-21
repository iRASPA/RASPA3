1) energy update during Monte-Carlo is minimal, just so to check energy-drift
   consists of: 
       framework-molecule VDW   
       framework-molecule Real  
       molecule-molecule VDW    
       molecule-molecule Real   
       Van der Waals (Tail)     
       Coulombic Ewald          
       intra VDW                
       intra Coulombic          
       polarization             
       dU/dlambda               
2) full energy composition and pressure computation are combined, sampling every 10 cycles

