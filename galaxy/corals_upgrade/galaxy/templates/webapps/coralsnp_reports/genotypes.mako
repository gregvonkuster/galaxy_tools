<%inherit file="/base.mako"/>
<%namespace file="/message.mako" import="render_msg" />

%if message:
    ${render_msg( message, 'done' )}
%endif

<div class="report">
    <div class="reportBody">
        <h3 align="center">
            All genotypes of uploaded samples
        </h3>
        <table align="center" class="colored">
            %if len(genotypes) == 0:
                <tr>
                    <td colspan="2">
                        There are no genotypes of uploaded samples
                    </td>
                </tr>
            %else:
                <tr class="header">
                    <td>Coral MLG Clonal ID</td>
                    <td>Coral MLG Rep Sample ID</td>
                    <td>Genetic Coral Species Call</td>
                    <td>Samples with this Genotype</td>
                </tr>
                <% ctr = 0 %>
                %for genotype in genotypes:
                    %if ctr % 2 == 1:
                        <tr class="odd_row">
                    %else:
                        <tr class="tr">
                    %endif
                        <td>${genotype[1]}</td>
                        <td>${genotype[2]}</td>
                        <td>${genotype[3]}</td>
                        <td><a href="${h.url_for( controller='samples', action='with_genotype', sort_id='default', order='default', genotype_id=genotype[0], coral_mlg_clonal_id=genotype[1], coral_mlg_rep_sample_id=genotype[2], genetic_coral_species_call=genotype[3] )}">Samples</a></td>
                    </tr>
                    <% ctr += 1 %>
                %endfor
            %endif
        </table>
    </div>
</div>

