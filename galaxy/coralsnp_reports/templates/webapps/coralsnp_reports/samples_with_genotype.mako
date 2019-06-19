<%inherit file="/base.mako"/>
<%namespace file="/message.mako" import="render_msg" />

%if message:
    ${render_msg( message, 'done' )}
%endif

<div class="report">
    <div class="reportBody">
        <h3 align="center">
            Samples with genotype having:
        </h3>
        <h4>
            coral mlg clonal id ${coral_mlg_clonal_id}<br/>
            coral mlg rep sample id&nbsp;${coral_mlg_rep_sample_id}<br/>
            genetic coral species call&nbsp${genetic_coral_species_call}<br/>
            bcoral genet id&nbsp${bcoral_genet_id}
        </h4>
        <table align="center" width="30%" class="colored">
            %if len(samples) == 0:
                <tr>
                    <td colspan="2">
                        There are no samples with this genotype
                    </td>
                </tr>
            %else:
                <tr class="header">
                    <td>affy id</td>
                    <td>sample id</td>
                    <td>genotype id</td>
                    <td>field call</td>
                    <td>collect date</td>
                    <td>user specimen id</td>
                    <td>registry id</td>
                    <td>depth</td>
                    <td>dna extract method</td>
                    <td>dna concentration</td>
                    <td>public after</td>
                    <td>% miss</td>
                    <td>% ref</td>
                    <td>% alt</td>
                    <td>% het</td>
                </tr>
                <% ctr = 0 %>
                %for sample in samples:
                    %if ctr % 2 == 1:
                        <tr class="odd_row">
                    %else:
                        <tr class="tr">
                    %endif
                        <td>${sample[0]}</td>
                        <td>${sample[1]}</td>
                        <td>${sample[2]}</td>
                        <td>${sample[3]}</td>
                        <td>${sample[4]}</td>
                        <td>${sample[5]}</td>
                        <td>${sample[6]}</td>
                        <td>${sample[7]}</td>
                        <td>${sample[8]}</td>
                        <td>${sample[9]}</td>
                        <td>${sample[10]}</td>
                        <td>${sample[11]}</td>
                        <td>${sample[12]}</td>
                        <td>${sample[13]}</td>
                        <td>${sample[14]}</td>
                    </tr>
                    <% ctr += 1 %>
                %endfor
            %endif
        </table>
    </div>
</div>
