<%inherit file="/base.mako"/>
<%namespace file="/message.mako" import="render_msg" />
<%namespace file="/sorting_base.mako" import="get_sort_url, get_css" />

%if message:
    ${render_msg( message, 'done' )}
%endif

${get_css()}

<div class="report">
    <div class="reportBody">
        <h3 align="center">Genotypes of samples uploaded per month</h3>
        <h4 align="center">
            Click the month to view the number of genotypes for samples uploaded per day of that month
        </h4>
        <table align="center" width="30%" class="colored">
            %if len(genotypes) == 0:
                <tr><td colspan="2">There are no genotypes for uploaded samples</td></tr>
            %else:
                <tr class="header">
                    <td class="half_width">
                        ${get_sort_url(sort_id, order, 'date', 'genotypes', 'per_month', 'Month')}
                        <span class='dir_arrow date'>${arrow}</span>
                    </td>
                    <td class="half_width">
                        ${get_sort_url(sort_id, order, 'num_genotypes', 'genotypes', 'per_month', 'Number of Genotypes')}
                        <span class='dir_arrow num_genotypes'>${arrow}</span>
                    </td>
                </tr>
                <% ctr = 0 %>
                %for sample in genotypes:
                    %if ctr % 2 == 1:
                        <tr class="odd_row">
                    %else:
                        <tr class="tr">
                    %endif
                        <td>
                            <a href="${h.url_for( controller='genotypes', action='specified_month', specified_date=sample[0]+'-01' )}">
                                ${sample[2]}&nbsp;${sample[3]}
                            </a>
                        </td>
                        <td>${sample[1]}</td>
                    </tr>
                    <% ctr += 1 %>
                %endfor
            %endif
        </table>
    </div>
</div>

