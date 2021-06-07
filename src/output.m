% get frame number and start preparing output
frame = floor(step/nop);
fprintf('\n*****  preparing output frame %d for %s \n',frame,runID);

% prepare for plotting
TX = {'Interpreter','Latex'}; FS = {'FontSize',18};
TL = {'TickLabelInterpreter','Latex'}; TS = {'FontSize',14};
UN = {'Units','Centimeters'};

axh = 6.00; axw = 7.50;
ahs = 1.00; avs = 1.00;
axb = 0.75; axt = 0.90;
axl = 1.75; axr = 0.90;

set(0,'DefaultFigureVisible',plot_op)

% prepare and plot figure for mechanical solution fields
fh1 = figure(1); clf; colormap(ocean);
fh = axb + 2*axh + 1*avs + axt;
fw = axl + 3*axw + 2*ahs + axr;
set(fh1,UN{:},'Position',[3 3 fw fh]);
set(fh1,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
set(fh1,'Color','w','InvertHardcopy','off');
set(fh1,'Resize','off');
ax(1) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+1*axh+1*avs axw axh]);
ax(2) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+1*axh+1*avs axw axh]);
ax(3) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+1*axh+1*avs axw axh]);
ax(4) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs axw axh]);
ax(5) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+0*axh+0*avs axw axh]);
ax(6) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+0*axh+0*avs axw axh]);

set(fh1, 'CurrentAxes', ax(1))
imagesc(xc,zc,-W(:      ,ic)-0.*WBG(:      ,ic)); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['matrix z-velocity $W$'],TX{:},FS{:}); set(gca,'XTickLabel',[]);
text(0,0.9*L,['time ',num2str(time,4)],TX{:},FS{:},'HorizontalAlignment','center','VerticalAlignment','middle')
set(fh1, 'CurrentAxes', ax(2))
imagesc(xc,zc, U(ic,:      )+0.*UBG(ic,:      )); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['matrix x-velocity $U$'],TX{:},FS{:}); set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
set(fh1, 'CurrentAxes', ax(3))
imagesc(xc,zc, P(ic,ic)); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['pore pressure $P$'],TX{:},FS{:}); set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
set(fh1, 'CurrentAxes', ax(4))
imagesc(xc,zc,-w(:      ,ic)); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['segr. z-velocity $w$'],TX{:},FS{:});
set(fh1, 'CurrentAxes', ax(5))
imagesc(xc,zc, u(ic,:      )); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['segr. x-velocity $u$'],TX{:},FS{:}); set(gca,'YTickLabel',[]);
set(fh1, 'CurrentAxes', ax(6))
imagesc(xc,zc, p(ic,ic)); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['comp. pressure $p$'],TX{:},FS{:}); set(gca,'YTickLabel',[]);
drawnow;

% prepare and plot figure for material properties
fh2 = figure(2); clf; colormap(ocean);
fh = axb + 2*axh + 1*avs + axt;
fw = axl + 3*axw + 2*ahs + axr;
set(fh2,UN{:},'Position',[6 6 fw fh]);
set(fh2,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
set(fh2,'Color','w','InvertHardcopy','off');
set(fh2,'Resize','off');
ax(1) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+1*axh+1*avs axw axh]);
ax(2) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+1*axh+1*avs axw axh]);
ax(3) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+1*axh+1*avs axw axh]);
ax(4) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs axw axh]);
ax(5) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+0*axh+0*avs axw axh]);
ax(6) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+0*axh+0*avs axw axh]);

set(fh2, 'CurrentAxes', ax(1))
imagesc(xc,zc,max(-6, log10(K(ic,ic)))); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['log$_{10}$ Darcy coefficient $K$'],TX{:},FS{:}); set(gca,'XTickLabel',[]);
text(0,0.9*L,['time ',num2str(time,4)],TX{:},FS{:},'HorizontalAlignment','center','VerticalAlignment','middle')
set(fh2, 'CurrentAxes', ax(2))
imagesc(xc,zc,      ups(ic,ic) ); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['decompaction rate $\dot{\upsilon}$'],TX{:},FS{:}); set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
set(fh2, 'CurrentAxes', ax(3))
imagesc(xc,zc,log10(eps(ic,ic))); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['log$_{10}$ strain rate $\dot{\varepsilon}$'],TX{:},FS{:}); set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
set(fh2, 'CurrentAxes', ax(4))
imagesc(xc,zc,log10( eta_vep(ic,ic))); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['log$_{10}$ shear viscosity $\eta$'],TX{:},FS{:});
set(fh2, 'CurrentAxes', ax(5))
imagesc(xc,zc,log10(zeta_vep(ic,ic))); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['log$_{10}$ plastic damage $\varepsilon_p$'],TX{:},FS{:}); set(gca,'YTickLabel',[]);
set(fh2, 'CurrentAxes', ax(6))
imagesc(xc,zc,log10(tau(ic,ic))); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['log$_{10}$ shear stress $\tau$'],TX{:},FS{:}); set(gca,'YTickLabel',[]);
drawnow;

% prepare and plot figure for thermo-chemical solution fields
fh3 = figure(3); clf; colormap(ocean);
fh = axb + 2*axh + 1*avs + axt;
fw = axl + 3*axw + 2*ahs + axr;
set(fh3,UN{:},'Position',[9 9 fw fh]);
set(fh3,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
set(fh3,'Color','w','InvertHardcopy','off');
set(fh3,'Resize','off');
ax(1) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+1*axh+1*avs axw axh]);
ax(2) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+1*axh+1*avs axw axh]);
ax(3) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+1*axh+1*avs axw axh]);
ax(4) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs axw axh]);
ax(5) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+0*axh+0*avs axw axh]);
ax(6) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+0*axh+0*avs axw axh]);

set(fh3, 'CurrentAxes', ax(1))
imagesc(xc,zc,T(ic,ic)-0.*mean(T(ic,ic),2)); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['Temperature $T$'],TX{:},FS{:}); set(gca,'XTickLabel',[]);
text(0,0.9*L,['time ',num2str(time,4)],TX{:},FS{:},'HorizontalAlignment','center','VerticalAlignment','middle')
set(fh3, 'CurrentAxes', ax(2))
imagesc(xc,zc,f(ic,ic)); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['Melt fraction $\phi$'],TX{:},FS{:}); set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
set(fh3, 'CurrentAxes', ax(3))
if Da>0
    imagesc(xc,zc,RctR_f(ic,ic)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['Melting Rate $\Gamma$'],TX{:},FS{:}); set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
else
    imagesc(xc,zc,fq(ic,ic)-f(ic,ic)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['Melt diseqilibrium $f^q-f$'],TX{:},FS{:}); set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
end
set(fh3, 'CurrentAxes', ax(4))
imagesc(xc,zc,MAJ(ic,ic)); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['Major elements $C_\mathrm{maj}$'],TX{:},FS{:});
set(fh3, 'CurrentAxes', ax(5))
imagesc(xc,zc,TRI(ic,ic)); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['Incomp. tracer $C_\mathrm{tri}$'],TX{:},FS{:}); set(gca,'YTickLabel',[]);
set(fh3, 'CurrentAxes', ax(6))
imagesc(xc,zc,TRC(ic,ic)); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['Comp. tracer $C_\mathrm{trc}$'],TX{:},FS{:}); set(gca,'YTickLabel',[]);
drawnow;


% prepare and plot figure for thermo-chemical solution fields
fh4 = figure(4); clf; colormap(ocean);
fh = axb + 2*axh + 1*avs + axt;
fw = axl + 2*axw + 1*ahs + axr;
set(fh4,UN{:},'Position',[12 12 fw fh]);
set(fh4,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
set(fh4,'Color','w','InvertHardcopy','off');
set(fh4,'Resize','off');
ax(1) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+1*axh+1*avs axw axh]);
ax(2) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+1*axh+1*avs axw axh]);
ax(3) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs axw axh]);
ax(4) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+0*axh+0*avs axw axh]);

set(fh4, 'CurrentAxes', ax(1))
imagesc(xc,zc,IRP(ic,ic)); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['Rad. parent isotope $C_\mathrm{irp}$'],TX{:},FS{:}); set(gca,'XTickLabel',[]);
text(0,0.9*L,['time ',num2str(time,4)],TX{:},FS{:},'HorizontalAlignment','center','VerticalAlignment','middle')
set(fh4, 'CurrentAxes', ax(2))
imagesc(xc,zc,IRD(ic,ic)); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['Rad. daughter isotope $C_\mathrm{ird}$'],TX{:},FS{:}); set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
set(fh4, 'CurrentAxes', ax(3))
imagesc(xc,zc,ISS(ic,ic)); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['Sld stable isotope $C_\mathrm{iss}$'],TX{:},FS{:});
set(fh4, 'CurrentAxes', ax(4))
imagesc(xc,zc,ISF(ic,ic)); axis ij equal tight; box on; cb = colorbar;
set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['Fld stable isotope $C_\mathrm{isf}$'],TX{:},FS{:}); set(gca,'YTickLabel',[]);
drawnow;


% prepare and plot figure for solution residuals
if plot_cv
    fh5 = figure(5); clf; colormap(ocean);
    fh = axb + 2*axh + 1*avs + axt;
    fw = axl + 3*axw + 2*ahs + axr;
    set(fh5,UN{:},'Position',[15 15 fw fh]);
    set(fh5,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
    set(fh5,'Color','w','InvertHardcopy','off');
    set(fh5,'Resize','off');
    ax(1) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+1*axh+1*avs axw axh]);
    ax(2) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+1*axh+1*avs axw axh]);
    ax(3) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+1*axh+1*avs axw axh]);
    ax(4) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs axw axh]);
    ax(5) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+0*axh+0*avs axw axh]);
    ax(6) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+0*axh+0*avs axw axh]);

    if bnchmrk
        set(fh5, 'CurrentAxes', ax(1))
        imagesc(xc,zc,-(W-W_mms)./(1e-16+norm(W_mms(:),2)./N)); axis ij equal tight; box on; cb = colorbar;
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['err. matrix z-velocity $W$'],TX{:},FS{:}); set(gca,'XTickLabel',[]); 
        set(fh5, 'CurrentAxes', ax(2))
        imagesc(xc,zc, (U-U_mms)./(1e-16+norm(U_mms(:),2)./N)); axis ij equal tight; box on; cb = colorbar;
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['err. matrix x-velocity $U$'],TX{:},FS{:}); set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
        set(fh5, 'CurrentAxes', ax(3))
        imagesc(xc,zc, (P-P_mms)./(1e-16+norm(P_mms(:),2)./N)); axis ij equal tight; box on; cb = colorbar;
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['err. pore pressure $P$'],TX{:},FS{:}); set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
        set(fh5, 'CurrentAxes', ax(4))
        imagesc(xc,zc, (T-T_mms)./(1e-16+norm(T_mms(:),2)./N)); axis ij equal tight; box on; cb = colorbar;
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['err. temperature $T$'],TX{:},FS{:});
        set(fh5, 'CurrentAxes', ax(5))
        imagesc(xc,zc, (C-C_mms)./(1e-16+norm(C_mms(:),2)./N)); axis ij equal tight; box on; cb = colorbar;
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['err. composition $C$'],TX{:},FS{:}); set(gca,'YTickLabel',[]);
        set(fh5, 'CurrentAxes', ax(6))
        imagesc(xc,zc, (f-f_mms)./(1e-16+norm(f_mms(:),2)./N)); axis ij equal tight; box on; cb = colorbar;
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['err. liquid fraction $\phi$'],TX{:},FS{:}); set(gca,'YTickLabel',[]);
        drawnow;
    else
        set(fh5, 'CurrentAxes', ax(1))
        imagesc(xc,zc,( FW./(1e-16+norm(W(:)+WBG(:),2)./N))); axis ij equal tight; box on; cb = colorbar;
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['res. matrix z-velocity $W$'],TX{:},FS{:}); set(gca,'XTickLabel',[]);
        set(fh5, 'CurrentAxes', ax(2))
        imagesc(xc,zc,(-FU./(1e-16+norm(W(:)+WBG(:),2)./N))); axis ij equal tight; box on; cb = colorbar;
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['res. matrix x-velocity $U$'],TX{:},FS{:}); set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
        set(fh5, 'CurrentAxes', ax(3))
        imagesc(xc,zc,(-FP./(1e-16+norm(P(:),2)       ./N))); axis ij equal tight; box on; cb = colorbar;
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['res. pore pressure $P$'],TX{:},FS{:}); set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]);
        set(fh5, 'CurrentAxes', ax(4))
        imagesc(xc,zc,(-res_T*(dt)./(1e-16+norm(T(:),2)       ./N))); axis ij equal tight; box on; cb = colorbar;
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['res. temperature $T$'],TX{:},FS{:});
        set(fh5, 'CurrentAxes', ax(5))
        imagesc(xc,zc,(-res_MAJ*(dt)./(1e-16+norm(MAJ(:),2)       ./N))); axis ij equal tight; box on; cb = colorbar;
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['res. composition $C$'],TX{:},FS{:}); set(gca,'YTickLabel',[]);
        set(fh5, 'CurrentAxes', ax(6))
        imagesc(xc,zc,(-res_f*(dt)./(1e-16+norm(f(:),2)       ./N))); axis ij equal tight; box on; cb = colorbar;
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['res. liquid fraction $\phi$'],TX{:},FS{:}); set(gca,'YTickLabel',[]);
        drawnow;
    end
    
    if ~bnchmrk
        figure(6); clf;
        pp = linspace(-1,max(Pe(:)),1e3);
        
        plot(pp,eps0.*ones(size(pp)),'k',pp,min(Ty+pp,2*Ty+pp/2),'r',Pe(:),yieldt(:),'r.','LineWidth',2); axis equal tight; box on; hold on;
        scatter(Pe(:),tau(:),20,(eta(:)),'filled'); colorbar; colormap(ocean);
        
        set(gca,'TickLabelInterpreter','latex','FontSize',15)
        title('Failure Criterion','Interpreter','latex','FontSize',22)
        xlabel('Effective Pressure','Interpreter','latex','FontSize',18)
        ylabel('Dev. Stress Magnitude','Interpreter','latex','FontSize',18)
        drawnow;
        
        fh7 = figure(7); clf;
        TT = linspace(0,1,1e3);
        [~,CCS,CCf] = equilibrium(TT,0.*TT,0.*TT,perT,perCs,perCf,clap,PhDg);
        
        plot(CCS,TT,'k-','LineWidth',2); axis tight; hold on; box on;
        plot(CCf,TT,'k-','LineWidth',2);
        
        plot([perCs,1],[0,0],'k-','LineWidth',1.5)
        plot([0,perCs],[perT,perT],'k-','LineWidth',1.5)
        plot([perCs,perCf],[perT,perT],'k-','LineWidth',1)
        plot([perCs,perCs],[-perT/2,perT],'k-','LineWidth',1.5)

        plot(MAJ,T-Pt*clap,'k.',MAJs,T-Pt*clap,'b.',MAJf,T-Pt*clap,'r.','LineWidth',2,'MarkerSize',20);

        set(gca,'XTick',[0:0.2:1],'YTick',[0:0.2:1],'TickLabelInterpreter','latex','FontSize',15)
        text(perCs/2,(perT-perT/2)/2,'sol 1 + sol 2','Interpreter','latex','FontSize',18,'HorizontalAlignment','center','VerticalAlignment','middle')
        text(perCs/2,perT+0.4*(1-perT),'sol 1 + liq','Interpreter','latex','FontSize',18,'HorizontalAlignment','center','VerticalAlignment','middle')
        text((perCs+1)/2,perT*0.4,'sol 2 + liq','Interpreter','latex','FontSize',18,'HorizontalAlignment','center','VerticalAlignment','middle')
        text((perCs+1)/2,-perT/4,'sol 2 + sol 3','Interpreter','latex','FontSize',18,'HorizontalAlignment','center','VerticalAlignment','middle')
        text(perCf-0.2,perT+0.75*(1-perT),'liq','Interpreter','latex','FontSize',18,'HorizontalAlignment','center','VerticalAlignment','middle')
        title('Phase Diagram','Interpreter','latex','FontSize',22)
        xlabel('Composition','Interpreter','latex','FontSize',18)
        ylabel('Temperature','Interpreter','latex','FontSize',18)
    end
end

% save output frame
if save_op
    % clean workspace
    clear Wi Ui Pi dWi dUi dPi Rfo RTo RMAJo RTRIo RTRCo RIRPo RIRDo RISSo RISFo Div_tz Div_tx dtW dtU dtP yield_GM yield_MC
    clear ax cb fw axw axh avs ahs axl axr axt axb fh fw TX TL FS TS UN p0 Vel k kk pp CCs CCf TT
    clear V_GrdT Div_fMAJV Div_fTRIV DIV_fTRCV Div_fIRPV Div_fIRDV V_GrdISS V_GrdISF V_GrdDMG
    clear Lpl_f Lpl_T Lpl_MAJ Lpl_TRI Lpl_TRC Lpl_IRP Lpl_IRD Lpl_ISS Lpl_ISF Lpl_DMG

    name = ['../out/',runID,'/',runID,'_svp_',num2str(frame)];
    print(fh1,name,'-dpng','-r300');
    name = ['../out/',runID,'/',runID,'_mat_',num2str(frame)];
    print(fh2,name,'-dpng','-r300');
    name = ['../out/',runID,'/',runID,'_stc_',num2str(frame)];
    print(fh3,name,'-dpng','-r300');
    name = ['../out/',runID,'/',runID,'_sis_',num2str(frame)];
    print(fh4,name,'-dpng','-r300');
    if ~bnchmrk
        name = ['../out/',runID,'/',runID,'_phs_',num2str(frame)];
        print(fh7,name,'-dpng','-r300');
    end
    
    clear fh1 fh2 fh3 fh4 fh5 fh7;

    name = ['../out/',runID,'/',runID,'_cont'];
    save([name,'.mat']);
    name = ['../out/',runID,'/',runID,'_',num2str(frame)];
    save([name,'.mat']);
    
    if step == 0
        logfile = ['../out/',runID,'/',runID,'.log'];
        if exist(logfile,'file'); delete(logfile); end
        diary(logfile)
    end
    
end
