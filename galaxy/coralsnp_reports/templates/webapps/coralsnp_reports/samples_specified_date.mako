<%inherit file="/base.mako"/>
<%namespace file="/message.mako" import="render_msg" />

%if message:
    ${render_msg( message, 'done' )}
%endif

<div class="report">
    <div class="reportBody">
        <h3 align="center">
            Samples for ${day_label},
            &nbsp;${month_label}&nbsp;${day_of_month},
            &nbsp;${year_label}
        </h3>
        <table align="center" width="30%" class="colored">
            %if len(samples) == 0:
                <tr>
                    <td colspan="2">
                        There are no samples for ${day_label},
                        &nbsp;${month_label}&nbsp;${day_of_month},
                        &nbsp;${year_label}
                    </td>
                </tr>
            %else:
                <tr class="header">
                    <td>affy id</td>
                    <td>sample id</td>
                    <td>genotype id</td>
                    <td>field call</td>
                    <td>collector last</td>
                    <td>first</td>
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
                        <td><a href="${h.url_for( controller='genotypes', action='for_sample', sort_id='default', order='default', genotype_id=sample[2], affy_id=sample[0], sample_id=sample[1], user_specimen_id=sample[7] )}">${sample[2]}</a></td>
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
                        <td>${sample[15]}</td>
                        <td>${sample[16]}</td>
                    </tr>
                    <% ctr += 1 %>
                %endfor
            %endif
        </table>
    </div>
</div>
