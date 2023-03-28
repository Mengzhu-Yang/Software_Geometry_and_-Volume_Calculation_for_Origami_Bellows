function [app.UIAxes2,app.UIAxes3,app.UIAxes4] = GeoandVolume(phi1,phi2,r,n,m,d1overr,value,Flat)
if strcmp(value, 'Miura')
                   text(app.UIAxes4,0,0,0.5,'There is no S3 for Miura');
                   if app.NumberofSolutions.Value == 0
                       text(app.UIAxes2,0.1,0.1,0.5,'There is no S1')
                       text(app.UIAxes3,0.1,0.1,0.5,'There is no S2') ;
                       %%%%%%%%%%
                   elseif app.NumberofSolutions.Value == 1
                       [Node, Panel, V1, Vs] = MiuraMonoS1(phi1,phi2,r,n,m,d1overr);
                       app.Vs.Value=Vs;
                       app.V1.Value=V1;
                       PlotOri(app.UIAxes2,Node,Panel,[],'PanelColor','g');
                       text(app.UIAxes3,0.1,0.1,0.5,'There is no S2') ;
                       axis(app.UIAxes2,'equal');
                       axis(app.UIAxes3,'equal');
                       axis(app.UIAxes4,'equal');

                       l2=axis(app.UIAxes2);
                       axis(app.UIAxes3,l2);
                       zlim(app.UIAxes3,l2(1,5:6));
                       axis(app.UIAxes4,l2);
                       zlim(app.UIAxes4,l2(1,5:6));
                       %%%%%%%%%%
                   elseif app.NumberofSolutions.Value == 2
                       if strcmp(Flat, 'Flat-foldable')
                           text(app.UIAxes2,0.1,0.1,0.5,'Flat-foldable')
                           app.Vs.Value=0;
                           app.V1.Value=0;
                       else
                           [Node, Panel, V1, Vs] = MiuraBiS1(phi1,phi2,r,n,m,d1overr);
                           app.Vs.Value=Vs;
                           app.V1.Value=V1;
                           PlotOri(app.UIAxes2,Node,Panel,[],'PanelColor','g');
                       end
                       [Node, Panel,V2,Vd] = MiuraBiS2(phi1,phi2,r,n,m,d1overr);
                       app.VdS2.Value=Vd;
                       app.V2.Value=V2;
                       PlotOri(app.UIAxes3,Node,Panel,[],'PanelColor','g');
                       axis(app.UIAxes2,'equal');
                       axis(app.UIAxes3,'equal');
                       axis(app.UIAxes4,'equal');

                       l3=axis(app.UIAxes3);
                       axis(app.UIAxes2,l3);
                       zlim(app.UIAxes2,l3(1,5:6));
                       axis(app.UIAxes4,l3);
                       zlim(app.UIAxes4,l3(1,5:6));
                   end


                   %%%%%Kresling
               elseif strcmp(value, 'Kresling')
                   if app.NumberofSolutions.Value == 0
                       text(app.UIAxes2,0.1,0.1,0.5,'There is no S1') ;
                       text(app.UIAxes3,0.1,0.1,0.5,'There is no S2') ;
                       text(app.UIAxes4,0.1,0.1,0.5,'There is no S3') ;
                       %%%%%%%%%%
                   elseif app.NumberofSolutions.Value == 2
                       [Node, Panel,V1,Vs] = KreslingMonoS1(phi1,phi2,r,n,m);
                       app.Vs.Value=Vs;
                       app.V1.Value=V1;
                       PlotOri(app.UIAxes2,Node,Panel,[],'PanelColor','[0.67, 0.83, 0.9]/256');

                       text(app.UIAxes3,0.1,0.1,0.5,'There is no S2') ;

                       [Node, Panel,V3,Vd] = KreslingS3(phi1,phi2,r,n,m);
                       app.VdS3.Value=Vd;
                       app.V3.Value=V3;
                       PlotOri(app.UIAxes4,Node,Panel,[],'PanelColor','[0.67, 0.83, 0.9]/256');
                       %%%%%%%%%%
                   elseif app.NumberofSolutions.Value == 3
                       if strcmp(Flat, 'Flat-foldable')
                           text(app.UIAxes2,0.1,0.1,0.5,'Flat-foldable')
                           app.Vs.Value=0;
                           app.V1.Value=0;
                       else
                           [Node, Panel,V1,Vs] = KreslingBiS1(phi1,phi2,r,n,m);
                           app.Vs.Value=Vs;
                           app.V1.Value=V1;
                           PlotOri(app.UIAxes2,Node,Panel,[],'PanelColor','[0.67, 0.83, 0.9]/256');
                       end
                       %
                       [Node, Panel,V2,Vd] = KreslingBiS2(phi1,phi2,r,n,m);
                       app.VdS2.Value=Vd;
                       app.V2.Value=V2;
                       PlotOri(app.UIAxes3,Node,Panel,[],'PanelColor','[0.67, 0.83, 0.9]/256');

                       [Node, Panel,V3,Vd] = KreslingS3(phi1,phi2,r,n,m);
                       app.VdS3.Value=Vd;
                       app.V3.Value=V3;
                       PlotOri(app.UIAxes4,Node,Panel,[],'PanelColor','[0.67, 0.83, 0.9]/256');
                   end
                   axis(app.UIAxes2,'equal');
                   axis(app.UIAxes3,'equal');
                   axis(app.UIAxes4,'equal');

                   l4=axis(app.UIAxes4);
                   axis(app.UIAxes2,l4);
                   zlim(app.UIAxes2,l4(1,5:6))
                   axis(app.UIAxes3,l4);
                   zlim(app.UIAxes3,l4(1,5:6))
               end
               axis(app.UIAxes2,'vis3d')
               axis(app.UIAxes3,'vis3d')
               axis(app.UIAxes4,'vis3d')
               
end