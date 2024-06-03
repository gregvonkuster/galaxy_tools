<%inherit file="/base.mako"/>
<%namespace file="/message.mako" import="render_msg" />

%if message:
    ${render_msg( message, 'done' )}
%endif

<div class="report">
    <div class="reportBody">
        <h3 align="center">
            All phenotypes of uploaded samples
        </h3>
        <table align="center" class="colored">
            %if len(phenotypes) == 0:
                <tr>
                    <td colspan="2">
                        There are no phenotypes of uploaded samples
                    </td>
                </tr>
            %else:
                <tr class="header">
                    <td>Disease Resist</td>
                    <td>Bleach Resist</td>
                    <td>Mortality</td>
                    <td>TLE</td>
                    <td>Spawning</td>
                    <td>Sperm Motility</td>
                    <td>Healing Time</td>
                    <td>Samples with this Phenotype</td>
                </tr>
                <% ctr = 0 %>
                %for phenotype in phenotypes:
                    %if ctr % 2 == 1:
                        <tr class="odd_row">
                    %else:
                        <tr class="tr">
                    %endif
                        <td>${phenotype[1]}</td>
                        <td>${phenotype[2]}</td>
                        <td>${phenotype[3]}</td>
                        <td>${phenotype[4]}</td>
                        <td>${phenotype[5]}</td>
                        <td>${phenotype[6]}</td>
                        <td>${phenotype[7]}</td>
                        <td><a href="${h.url_for( controller='samples', action='with_phenotype', sort_id='default', order='default', phenotype_id=phenotype[0], disease_resist=phenotype[1], bleach_resist=phenotype[2], mortality=phenotype[3], tls=phenotype[4], spawning=phenotype[5], sperm_motility=phenotype[6], healing_time=phenotype[7] )}">Samples</a></td>
                    </tr>
                    <% ctr += 1 %>
                %endfor
            %endif
        </table>
    </div>
</div>

